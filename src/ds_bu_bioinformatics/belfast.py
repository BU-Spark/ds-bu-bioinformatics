"""Belfast / FASTMAN GSE116918: Almac sense probes with metastasis and BCR outcomes after radiotherapy."""
import csv
import io
import re

import numpy as np
import pandas as pd

from .cambridge import read_matrix
from .expression import open_text
from .geo_cohort_pipeline import add_outcomes, decode_bcr

SERIES = 'GSE116918'
PLATFORM = 'GPL25318'
STUDY = 'Belfast / FASTMAN'
SNAPSHOT = 'belfast-v1'
SCALE = 'Deposited RMA probeset-level log2 intensity (GEO data processing); no added transformation'
MATRIX_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116918/matrix/GSE116918_series_matrix.txt.gz'
PLATFORM_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL25318&targ=self&form=text&view=full'

PRIMARY_ENDPOINT = 'metastasis'
# endpoint -> (GEO event field, GEO time field, definition)
ENDPOINTS = {
    'metastasis': ('geo__met_event_1_yes_0_no', 'geo__follow_up_time_met_months',
                   'Metastasis after radical radiotherapy with ADT; GEO "met event (1=yes, 0=no)" and '
                   '"follow-up time (met, months)"'),
    'bcr': ('geo__bcr_event_1_yes_0_no', 'geo__follow_up_time_bcr_months',
            'Biochemical recurrence after radical radiotherapy with ADT (not post-prostatectomy BCR); GEO '
            '"bcr event (1=yes, 0=no)" and "follow-up time (bcr, months)"'),
}
CODES = {'1': 1, '0': 0}
TISSUE = 'Primary Prostate Cancer'
TITLE = re.compile(r'S\d+-\d+')

SENSE = {'Sense (Fully Exonic)', 'Sense (Includes Intronic)'}
NO_MATCH = {'No Genome Match', 'No Transcript Match'}
ORIENTATIONS = SENSE | NO_MATCH | {'AntiSense'}
PLATFORM_COLUMNS = ['ID', 'Orientation', 'Gene Symbol', 'Entrez Gene', 'Ensembl Gene ID']


def endpoint_outcome(name):
    event_field, time_field, _ = ENDPOINTS[name]
    def outcome(record):
        return decode_bcr(record.get(event_field), record.get(time_field), CODES)
    return outcome


outcome = endpoint_outcome('bcr')   # shared bcr_* columns, validated by geo_cohort_pipeline


def add_endpoint(clinical, name):
    """Add <name>_event/_months/_eligibility_reason/_eligible/_time_unit using the shared decoding rules."""
    decoded = add_outcomes(clinical, endpoint_outcome(name))
    if name == 'bcr':
        return decoded
    renamed = {f'bcr_{field}': f'{name}_{field}' for field in ('event', 'months', 'eligibility_reason', 'eligible', 'time_unit')}
    return clinical.join(decoded[list(renamed)].rename(columns=renamed))


def clinical_tables(metadata):
    audit = metadata.copy()
    required = {'geo__tissue', *(field for e in ENDPOINTS.values() for field in e[:2])}
    if not required.issubset(audit.columns):
        raise ValueError(f'Belfast GEO characteristics missing: {sorted(required - set(audit.columns))}')
    primary = audit.geo__tissue.eq(TISSUE) & audit.title.str.fullmatch(TITLE)
    if not primary.all():
        raise ValueError('Unrecognized Belfast tissue/title; every deposited profile should be a primary biopsy')
    audit['patient_id'] = audit.title
    audit['selected'] = primary
    audit['selection_reason'] = np.where(primary, 'primary_tumor', 'not_primary_tumor')
    clinical = audit.loc[primary].reset_index(drop=True)
    if clinical.empty or not clinical.patient_id.is_unique:
        raise ValueError('Empty/duplicate primary Belfast subjects')
    for name in ('bcr', PRIMARY_ENDPOINT):
        clinical = add_endpoint(clinical, name)
    clinical['outcome_clinical_source'] = 'GEO GSE116918 sample characteristics'
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title']].copy()
    return clinical, samples, audit


def read_platform(path):
    with open_text(path) as handle:
        text = handle.read()
    if f'^PLATFORM = {PLATFORM}' not in text or '!platform_table_begin' not in text:
        raise ValueError(f'{PLATFORM} table absent')
    table = text.split('!platform_table_begin', 1)[1].lstrip('\r\n')
    if '!platform_table_end' not in table:
        raise ValueError(f'Truncated {PLATFORM} annotation')
    table = table.split('!platform_table_end', 1)[0]
    frame = pd.read_csv(io.StringIO(table), sep='\t', dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
    if not set(PLATFORM_COLUMNS).issubset(frame.columns):
        raise ValueError(f'Unexpected {PLATFORM} schema')
    return frame[PLATFORM_COLUMNS]


def platform_annotations(path, expression):
    """One row per expression probe. Only sense-strand probes with one Entrez ID and one symbol are eligible."""
    frame = read_platform(path)
    if frame.ID.eq('').any() or not frame.ID.is_unique:
        raise ValueError(f'Missing/duplicate {PLATFORM} probe IDs')
    if not set(expression.index).issubset(frame.ID):
        raise ValueError(f'Expression probes missing from {PLATFORM}')
    unknown = set(frame.Orientation) - ORIENTATIONS
    if unknown:
        raise ValueError(f'Unrecognized {PLATFORM} orientation(s): {sorted(unknown)}')
    a = frame.set_index('ID').loc[expression.index].reset_index()
    a.columns = ['feature_id', 'orientation', 'gene_symbol', 'gene_id', 'ensembl_gene_id']
    sense = a.orientation.isin(SENSE)
    blank = a.gene_symbol.isin(['', '---']) | a.gene_id.isin(['', '---'])
    multi = a.gene_symbol.str.contains('///') | a.gene_id.str.contains('///')
    single = sense & ~blank & ~multi & a.gene_id.str.fullmatch(r'\d+') & a.gene_symbol.str.fullmatch(r'[^;/,\s]+')
    pairs = a.loc[single, ['gene_id', 'gene_symbol']].drop_duplicates()
    bad_ids = set(pairs.loc[pairs.gene_id.duplicated(False), 'gene_id'])
    bad_symbols = set(pairs.loc[pairs.gene_symbol.duplicated(False), 'gene_symbol'])
    conflict = single & (a.gene_id.isin(bad_ids) | a.gene_symbol.isin(bad_symbols))
    a['mapping_status'] = np.select(
        [a.orientation.eq('AntiSense'), a.orientation.isin(NO_MATCH), blank, multi, conflict, single],
        ['antisense_probe', 'no_transcript_match', 'unannotated', 'multi_gene', 'gene_symbol_identity_conflict', 'mapped'],
        default='invalid_annotation')
    a['annotation_eligible'] = a.mapping_status.eq('mapped')
    return a[['feature_id', 'gene_id', 'gene_symbol', 'ensembl_gene_id', 'orientation', 'mapping_status', 'annotation_eligible']]


def prepare(cache):
    meta, expression, processing = read_matrix(cache.fetch(f'{SERIES}_series_matrix.txt.gz', MATRIX_URL),
                                               series=SERIES, platform=PLATFORM)
    clinical, samples, audit = clinical_tables(meta)
    source = expression.loc[:, samples.sample_id]
    annotations = platform_annotations(cache.fetch(f'{PLATFORM}.soft.txt', PLATFORM_URL), source)
    details = {
        'sample_selection': f'GEO tissue "{TISSUE}" with S####-##### title; every deposited profile is one diagnostic '
                            'biopsy from one patient, selected independently of outcomes',
        'patient_identity': 'Patient label is the GEO sample title; no external registry linkage',
        'annotation': f'Frozen {PLATFORM} (Almac ADXPCv1a, Ensembl release 75): only Sense (Fully Exonic) or Sense '
                      '(Includes Intronic) probes with one Entrez ID and one symbol are mapped. AntiSense probes are '
                      'excluded even when the annotation names a gene, because they do not measure that gene\'s '
                      'sense transcript. Multi-gene, unannotated, no-match and identity-conflict probes are excluded.',
        'outcomes': {name: definition for name, (_, _, definition) in ENDPOINTS.items()},
        'primary_endpoint': PRIMARY_ENDPOINT,
        'definition_source': 'https://pubmed.ncbi.nlm.nih.gov/29045551/ (Jain et al. 2018)',
        'limitations': ['Diagnostic FFPE biopsy, not prostatectomy tissue', 'Treatment is radical radiotherapy with ADT',
                        'T stage is clinical; Gleason is a total score without primary/secondary patterns',
                        'GEO event counts (BCR 56, metastasis 22) differ from Jain 2018 (65, 24); deposited values are used',
                        'No raw CEL quality control', 'Deposited platform annotation (Ensembl 75) may be outdated']}
    return clinical, samples, audit, source, annotations, {}, processing, details


def summarize(clinical, audit, mapping, extra):
    return {}
