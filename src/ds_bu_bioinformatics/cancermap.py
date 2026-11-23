"""CancerMap GSE94767: human malignant regions, one explicit source donor sample."""
import numpy as np
import pandas as pd

from .exon_arrays import source_annotations
from .geo_cohort_pipeline import add_outcomes
from .geo_io import read_geo

SERIES = 'GSE94767'
PLATFORM = 'GPL5175'
STUDY = 'CancerMap'
SNAPSHOT = 'cancermap-v1'
SCALE = 'Deposited normalized gene-target expression; no added transformation or aggregation'
MATRIX_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94767/matrix/GSE94767_series_matrix.txt.gz'


def outcome(record):
    return None, None, 'missing_status'


ENDPOINTS = {'bcr': outcome}


def source_column(frame, names, label):
    present = [c for c in names if c in frame]
    if not present:
        raise ValueError(f'CancerMap {label} field absent; source columns are {list(frame)}')
    value = frame[present[0]].fillna('').astype(str).str.strip()
    for column in present[1:]:
        other = frame[column].fillna('').astype(str).str.strip()
        conflicts = value.ne('') & other.ne('') & value.ne(other)
        if conflicts.any():
            raise ValueError(f'Conflicting CancerMap {label} fields')
        value = value.mask(value.eq(''), other)
    return value


def clinical_tables(metadata):
    audit = metadata.copy()
    # Deposited tissue is the organ (Prostate), not the malignancy class.
    material = source_column(audit, ['geo__material_type'], 'material_type')
    donor = source_column(audit, ['geo__donor_id', 'geo__donor', 'geo__patient_id', 'geo__patient', 'geo__individual'], 'donor')
    human = audit.get('organism', pd.Series('', index=audit.index)).eq('Homo sapiens')
    malignant = material.str.lower().isin({'malignant', 'malignant tissue', 'tumour', 'tumor', 'cancer',
        'prostate cancer', 'prostate tumour', 'prostate tumor', 'tumour tissue', 'tumor tissue'})
    benign = material.str.lower().isin({'benign', 'benign tissue', 'normal', 'normal tissue', 'morphologically benign tissue',
                                     'prostate benign', 'benign prostate tissue'})
    source_name = audit.get('source_name', pd.Series('', index=audit.index)).fillna('').astype(str)
    xenograft = ~human | source_name.str.contains(r'xenograft|\bPDX\b', case=False, regex=True)
    if audit.empty or not audit.platform_id.eq(PLATFORM).all() or not (malignant | benign | xenograft).all():
        raise ValueError('Unrecognized CancerMap specimen; review source tissue/organism values')
    primary = human & malignant & ~xenograft
    if not primary.any() or donor[primary].isin(['', 'NA', 'N/A', 'UNKNOWN']).any():
        raise ValueError('CancerMap primary donor identity missing')
    audit['patient_id'] = SERIES + ':' + donor
    # No independent per-region quality ranking is supplied: a deterministic,
    # outcome-independent rule retains one original region without averaging.
    winners = set(audit.loc[primary].sort_values('sample_id').drop_duplicates('patient_id').sample_id)
    audit['selected'] = audit.sample_id.isin(winners)
    audit['selection_reason'] = np.select([audit.selected, primary, benign & human & ~xenograft, xenograft],
        ['primary_region_lowest_gsm', 'other_primary_region_same_donor', 'benign', 'nonhuman_or_xenograft'], default='unrecognized')
    clinical = add_outcomes(audit.loc[audit.selected].reset_index(drop=True), outcome)
    clinical['bcr_clinical_source'] = 'BCR status/duration not verified from deposited CancerMap metadata'
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title']].copy()
    return clinical, samples, audit


def prepare(cache):
    meta, expression, processing = read_geo(cache.fetch(SERIES + '_series_matrix.txt.gz', MATRIX_URL),
                                           series=SERIES, platforms={PLATFORM})
    clinical, samples, audit = clinical_tables(meta)
    source = expression.loc[:, samples.sample_id]
    genes, candidates = source_annotations(cache, PLATFORM, source)
    details = {'sample_selection': 'Human malignant tissue; one donor region chosen by lowest GSM identifier, independently of outcomes; all regions audited',
        'patient_identity': 'Explicit source donor characteristic; never infer one patient per assay',
        'annotation': 'Single-symbol source gene_assignment anchored to cached human RefSeq mapping; deposited target rows preserved',
        'expression_kind': 'log_microarray',
        'endpoints': {'bcr': {'definition': 'Biochemical recurrence', 'duration_available': False, 'unit': 'months'}},
        'study_source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94767',
        'limitations': ['Full-source schema/count verification pending; synthetic tests only',
            'Lowest GSM is a reproducible sampling rule, not an assertion that the selected region is clinically dominant; assess region sensitivity',
            'Deposited transcript-level summaries do not recover raw individual probes',
            'External patient-linked clinical supplement required for survival; missing outcomes do not remove expression subjects']}
    return clinical, samples, audit, source, genes, {'gene_mapping_candidates.csv.gz': candidates}, processing, details


def summarize(clinical, audit, mapping, extra):
    return {'survival_ready': bool(clinical.bcr_eligible.any()),
            'other_primary_regions_excluded': int(audit.selection_reason.eq('other_primary_region_same_donor').sum())}
