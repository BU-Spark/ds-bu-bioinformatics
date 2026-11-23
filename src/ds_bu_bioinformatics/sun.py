"""Sun/Goodison GSE25136: expression and recurrence status, without invented times."""
from .cambridge import read_matrix
from .geo_cohort_pipeline import add_outcomes, decode_bcr
from .geo_platform import annotations

SERIES = 'GSE25136'
PLATFORM = 'GPL96'
STUDY = 'Sun / Goodison'
SNAPSHOT = 'sun-v1'
SCALE = 'Deposited MAS5 globally scaled microarray signal; not assumed log2; no added transformation'
MATRIX_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25136/matrix/GSE25136_series_matrix.txt.gz'
PLATFORM_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96&targ=self&form=text&view=full'


def outcome(record):
    # GEO exposes a class label. Neither cohort follow-up nor literature medians
    # can supply individual censoring times.
    return decode_bcr(record.get('geo__recurrence_status'), None,
                      {'Recurrent': 1, 'Non-Recurrent': 0})


ENDPOINTS = {'recurrence': outcome}


def clinical_tables(metadata):
    required = {'geo__tissue', 'geo__recurrence_status', 'title', 'sample_id', 'geo_accession', 'platform_id'}
    if not required.issubset(metadata):
        raise ValueError('Sun: missing source tissue, identity or recurrence column')
    audit = metadata.copy()
    if (audit.empty or not audit.geo__tissue.eq('Prostate cancer primary tumor').all()
            or not audit.title.str.fullmatch(r'Prostate cancer primary tumor PR\d+').all()
            or not audit.platform_id.eq(PLATFORM).all()):
        raise ValueError('Sun: unrecognized specimen/platform/title; review selection')
    audit['patient_id'] = audit.title.str.rsplit(' ', n=1).str[-1]
    if not audit.patient_id.is_unique or not audit.sample_id.is_unique:
        raise ValueError('Sun: duplicate patient or assay')
    audit['selected'] = True
    audit['selection_reason'] = 'primary_tumor'
    clinical = add_outcomes(audit.reset_index(drop=True), outcome, 'recurrence')
    clinical['recurrence_clinical_source'] = 'GEO GSE25136 recurrence status only; no verified patient-level duration'
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title']].copy()
    return clinical, samples, audit


def prepare(cache):
    meta, expression, processing = read_matrix(cache.fetch('GSE25136_series_matrix.txt.gz', MATRIX_URL),
                                               series=SERIES, platform=PLATFORM)
    clinical, samples, audit = clinical_tables(meta)
    source = expression.loc[:, samples.sample_id]
    genes = annotations(cache.fetch('GPL96.soft.txt', PLATFORM_URL), source,
        platform=PLATFORM, gene_column='ENTREZ_GENE_ID', gene_pattern=r'\d+',
        species_column='Species Scientific Name', extra_columns={'RefSeq Transcript ID': 'refseq_original'})
    details = {
        'expression_kind': 'linear_microarray',
        'sample_selection': 'Require primary tumor tissue and unique PR* title-derived subjects; independent of recurrence class',
        'patient_identity': 'GEO-title-derived PR label; no independent registry linkage',
        'annotation': 'Single human Entrez/symbol from frozen GPL96; ambiguous/conflicting mappings quarantined; separate probes retained',
        'endpoints': {'recurrence': {'definition': 'Deposited recurrent/non-recurrent class; not relabeled BCR or DFS',
            'source_fields': ['recurrence status'], 'unit': 'months', 'duration_available': False,
            'missing_policy': 'All times nullable; all records ineligible for survival; never infer time from class'}},
        'study_source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25136',
        'limitations': ['Full-source acquisition/count verification pending; tested with synthetic source-shaped inputs',
            'No verified individual recurrence/censoring times: Cox, Kaplan-Meier and C-index unavailable',
            'Recurrence labels support class comparisons, not prospective risk calibration',
            'MAS5 values retained; expression normalization must be reviewed for downstream module comparisons',
            'Frozen platform labels may be outdated; no raw CEL reprocessing']}
    return clinical, samples, audit, source, genes, {}, processing, details


def summarize(clinical, audit, mapping, extra):
    ready = bool(clinical.recurrence_eligible.any())
    return {'survival_ready': ready, 'survival_blocker': None if ready else 'No verified patient-linked recurrence or censoring durations'}
