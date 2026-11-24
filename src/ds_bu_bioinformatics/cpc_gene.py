"""CPC-GENE GSE107299 processed Entrez matrix with exact source-patient linkage."""
import re
import csv

import numpy as np
import pandas as pd

from .gene_identifiers import annotate
from .geo_cohort_pipeline import add_outcomes
from .geo_io import align_columns, read_geo, read_tsv
from .expression import open_text

SERIES = 'GSE107299'
PANELS = ('GPL16686', 'GPL17586')
PLATFORM = ';'.join(PANELS)
ACCESSIONS = [SERIES, *PANELS]
STUDY = 'CPC-GENE'
SNAPSHOT = 'cpc-gene-v1'
SCALE = 'Deposited RMA/Brainarray Entrez summaries with source batch correction; no new transformation/aggregation'
BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107299'
MATRIX_URLS = {p: BASE + f'/matrix/{SERIES}-{p}_series_matrix.txt.gz' for p in PANELS}
EXPRESSION_URL = BASE + '/suppl/GSE107299_Matrix_processed_data.tsv.gz'


def outcome(record):
    return None, None, 'missing_status'


ENDPOINTS = {'bcr': outcome}


def patient_label(title):
    match = re.fullmatch(r'prostate tumour (CPCG\d+(?:_rep)?)', title)
    if not match:
        raise ValueError('Unrecognized CPC-GENE source patient title')
    return match.group(1)


def processed_matrix(path):
    """The deposited TSV has seven gene-annotation columns before assay values."""
    with open_text(path) as handle:
        header = next(csv.reader(handle, delimiter='\t'), [])
    annotation_fields = ['GeneID', 'Symbol_UCSC', 'Name_UCSC', 'Chr_UCSC',
                         'Start_UCSC', 'End_UCSC', 'RefSeq_UCSC']
    annotated = header[:7] == annotation_fields
    n = 7 if annotated else 1
    if (len(header) <= n or len(header) != len(set(header))
            or header[0] not in {'GeneID', 'Entrez'}
            or not all(re.fullmatch(r'CPCG\d+|GSM\d+', c) for c in header[n:])):
        raise ValueError('Unexpected CPC-GENE processed matrix schema')
    frame = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    if frame.empty or frame.iloc[:, 0].eq('').any() or frame.iloc[:, 0].duplicated().any():
        raise ValueError('CPC-GENE processed gene identifiers missing or duplicated')
    annotations = frame.iloc[:, :n].copy().rename(columns={header[0]: 'feature_id'})
    source = frame.iloc[:, n:].apply(lambda s: pd.to_numeric(s.mask(s.isin(['', 'NA', 'N/A', 'NaN'])), errors='raise'))
    if np.isinf(source.to_numpy()).any():
        raise ValueError('Infinite CPC-GENE source expression')
    source.index = pd.Index(frame.iloc[:, 0], name='feature_id')
    return source, annotations


def prepare(cache):
    metadata, processing = [], []
    for platform in PANELS:
        meta, _, methods = read_geo(cache.fetch(f'{SERIES}-{platform}_series_matrix.txt.gz', MATRIX_URLS[platform]),
                                    series=SERIES, platforms={platform}, metadata_only=True)
        metadata.append(meta)
        processing.extend(methods)
    meta = pd.concat(metadata, ignore_index=True)
    if not meta.sample_id.is_unique or not meta.geo__tissue.eq('prostate tumour').all():
        raise ValueError('CPC-GENE duplicate assay or unexpected tissue')
    meta['source_expression_id'] = meta.title.map(patient_label)
    meta['patient_id'] = meta.source_expression_id.str.replace(r'_rep$', '', regex=True)
    source, deposited_annotations = processed_matrix(cache.fetch('GSE107299_Matrix_processed_data.tsv.gz', EXPRESSION_URL))
    source, linkage = align_columns(source, meta, lambda r: [r['sample_id'], r['source_expression_id']])
    present = set(source.columns)
    # The GEO record explicitly states that CPCG0398_rep was averaged into
    # CPCG0398 upstream. Never average it a second time or count it as a patient.
    is_rep = meta.source_expression_id.str.endswith('_rep')
    if (is_rep & meta.source_expression_id.ne('CPCG0398_rep')).any():
        raise ValueError('Unexpected CPC-GENE replicate; review source processing')
    if any(x not in present for x in meta.loc[~is_rep, 'sample_id']):
        raise ValueError('Processed CPC-GENE matrix omits a non-replicate GEO subject')
    if any(x in present for x in meta.loc[is_rep, 'sample_id']):
        raise ValueError('Unexpected separate technical replicate in source processed matrix')
    meta['selected'] = meta.sample_id.isin(present)
    meta['selection_reason'] = np.where(meta.selected, 'deposited_processed_patient', 'technical_replicate_already_averaged_by_source')
    clinical = add_outcomes(meta.loc[meta.selected].reset_index(drop=True), outcome)
    if clinical.empty or not clinical.patient_id.is_unique:
        raise ValueError('CPC-GENE processed matrix contains duplicate patients')
    clinical['bcr_clinical_source'] = 'No verified patient-linked BCR table joined; source clinical covariates retained'
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title', 'source_expression_id']].copy()
    source = source.loc[:, samples.sample_id]
    genes, candidates = annotate(cache, source.index.tolist(), scope='entrezgene')
    details = {'sample_selection': 'Exact CPCG/GSM keys link processed columns to GEO patients; no positional join or new replicate average',
        'patient_identity': 'Source CPCG title label; preserve source batch/platform',
        'annotation': 'Brainarray Entrez source IDs validated through cached exact human Entrez mapping; original row IDs retained',
        'expression_kind': 'log_microarray',
        'endpoints': {'bcr': {'definition': 'Biochemical recurrence', 'duration_available': False, 'unit': 'months'}},
        'study_source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107299',
        'limitations': ['Full-source schema/count verification pending; synthetic tests only',
            'Source already summarizes probes, corrects batches and averages CPCG0398 technical replicates; these operations are not undone',
            'The clinical supplement/proteomic subset must be matched by source CPCG IDs before survival analysis',
            'Do not assume 55/76 proteomic-study cases equal all deposited expression patients']}
    return clinical, samples, meta, source, genes, {'expression_linkage.csv': linkage,
        'deposited_gene_annotations.csv': deposited_annotations,
        'gene_mapping_candidates.csv.gz': candidates}, sorted(set(processing)), details


def summarize(clinical, audit, mapping, extra):
    return {'survival_ready': bool(clinical.bcr_eligible.any()), 'source_averaged_replicate': 'CPCG0398_rep'}
