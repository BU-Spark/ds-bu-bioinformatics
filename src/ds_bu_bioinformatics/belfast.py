"""Belfast GSE116918 diagnostic biopsies: deposited RMA, BCR and metastasis."""
from .cambridge import read_matrix
from .geo_cohort_pipeline import add_outcomes, decode_bcr
from .geo_platform import annotations

SERIES = 'GSE116918'
PLATFORM = 'GPL25318'
STUDY = 'Belfast / FASTMAN (diagnostic biopsy; radiotherapy + ADT)'
SNAPSHOT = 'belfast-v1'
SCALE = 'Deposited probeset-level RMA-normalized signal; no additional transformation'
MATRIX_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116918/matrix/GSE116918_series_matrix.txt.gz'
PLATFORM_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL25318&targ=self&form=text&view=full'


def outcome(record):
    return decode_bcr(record.get('geo__bcr_event_1_yes_0_no'),
                      record.get('geo__follow_up_time_bcr_months'), {'1': 1, '0': 0})


def metastasis_outcome(record):
    return decode_bcr(record.get('geo__met_event_1_yes_0_no'),
                      record.get('geo__follow_up_time_met_months'), {'1': 1, '0': 0})


ENDPOINTS = {'bcr': outcome, 'metastasis': metastasis_outcome}


def clinical_tables(metadata):
    required = {'geo__tissue', 'title', 'sample_id', 'geo_accession', 'platform_id',
                'geo__bcr_event_1_yes_0_no', 'geo__follow_up_time_bcr_months',
                'geo__met_event_1_yes_0_no', 'geo__follow_up_time_met_months'}
    if not required.issubset(metadata):
        raise ValueError('Belfast: missing tissue, identity or endpoint source columns')
    audit = metadata.copy()
    if (audit.empty or not audit.geo__tissue.eq('Primary Prostate Cancer').all()
            or not audit.title.str.fullmatch(r'S1236-\d+').all()
            or not audit.platform_id.eq(PLATFORM).all()):
        raise ValueError('Belfast: unrecognized specimen, platform or subject title; review selection')
    audit['patient_id'] = audit.title
    if not audit.patient_id.is_unique or not audit.sample_id.is_unique:
        raise ValueError('Belfast: duplicate patient or GEO assay; no arbitrary replicate selection')
    audit['selected'] = True
    audit['selection_reason'] = 'primary_diagnostic_biopsy'
    clinical = audit.reset_index(drop=True)
    for prefix, decode in ENDPOINTS.items():
        clinical = add_outcomes(clinical, decode, prefix)
        clinical[f'{prefix}_clinical_source'] = f'GEO {SERIES} endpoint-specific sample event / follow-up time (months)'
    clinical['treatment_context'] = 'Cohort treated with radical radiotherapy and ADT; not patient-specific treatment adjudication'
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title']].copy()
    return clinical, samples, audit


def prepare(cache):
    meta, expression, processing = read_matrix(cache.fetch('GSE116918_series_matrix.txt.gz', MATRIX_URL),
                                               series=SERIES, platform=PLATFORM)
    clinical, samples, audit = clinical_tables(meta)
    source = expression.loc[:, samples.sample_id]
    genes = annotations(cache.fetch('GPL25318.soft.txt', PLATFORM_URL), source,
        platform=PLATFORM, gene_column='Ensembl Gene ID', gene_pattern=r'ENSG\d+',
        extra_columns={'Entrez Gene': 'entrez_original'})
    details = {
        'expression_kind': 'log_microarray',
        'sample_selection': 'Require primary-prostate-cancer tissue and unique S1236-* source subject titles; retain all regardless of clinical eligibility',
        'patient_identity': 'GEO title-derived subject label; no independent registry linkage',
        'annotation': 'Single deposited Ensembl gene ID and symbol; ambiguous/conflicting mappings quarantined; distinct probes never averaged',
        'endpoints': {
            'bcr': {'definition': 'PSA rise greater than 2 ng/mL above nadir followed by a further rise',
                    'time_origin': 'Start of external beam radiotherapy', 'censoring': 'Last follow-up',
                    'deposited_duration_contract_verified': False, 'unit': 'months',
                    'source_fields': ['bcr event (1=yes, 0=no)', 'follow-up time (bcr, months)']},
            'metastasis': {'definition': 'Radiologically documented metastasis',
                    'time_origin': 'Start of external beam radiotherapy', 'censoring': 'Last follow-up',
                    'deposited_duration_contract_verified': False, 'unit': 'months',
                    'source_fields': ['met event (1=yes, 0=no)', 'follow-up time (met, months)']}},
        'study_source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116918',
        'limitations': ['Full-source acquisition/count verification pending; tested with synthetic source-shaped inputs',
            'Diagnostic biopsies followed by RT/ADT differ from radical-prostatectomy cohorts; do not pool endpoints blindly',
            'GEO and publication event counts may differ; compute from cached records, never force literature counts',
            'Endpoint-specific deposited durations used as supplied; no independent adjudication of clinical time origin',
            'Frozen platform gene labels may be outdated; no raw CEL QC or reprocessing']}
    return clinical, samples, audit, source, genes, {}, processing, details


def summarize(clinical, audit, mapping, extra):
    return {}
