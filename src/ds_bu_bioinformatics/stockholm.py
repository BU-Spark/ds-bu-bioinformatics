"""Stockholm GSE70769 specimen rules; reuse shared Illumina/BCR validation."""
from . import cambridge

SERIES = 'GSE70769'
STUDY = 'Stockholm'
SNAPSHOT = 'stockholm-v1'
CLINICAL_COHORT = 'stockholm'
SAMPLE_POLICY = 'Require primary RP tumor title with STKHLM subject ID and publication tumor/RP identity; retain missing clinical outcomes'
MATRIX_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70769/matrix/GSE70769_series_matrix.txt.gz'


def read_matrix(path):
    return cambridge.read_matrix(path, series=SERIES)


def clinical_tables(metadata, publication):
    audit = metadata.copy()
    # GEO contains both the correctly spelled form and its deposited typo.
    prefixes = ('tumour tissue_robotic radical prostatectomy_',
                'tumour tissue_robotic radical prostatetctomy_')
    if not audit.title.str.startswith(prefixes).all():
        raise ValueError('Unrecognized Stockholm specimen title')
    audit['patient_id'] = audit.title.str.rsplit('_', n=1).str[-1]
    if not audit.patient_id.str.fullmatch(r'STKHLM\d+').all() or not audit.patient_id.is_unique:
        raise ValueError('Invalid/duplicate Stockholm subject IDs')
    if 'geo__sample_type' in audit and not audit.geo__sample_type.eq('Tumour').all():
        raise ValueError('Unexpected Stockholm GEO sample type')
    matched = publication.set_index('patient_id').reindex(audit.patient_id)
    if matched.paper__subject_id.isna().any():
        raise ValueError('Stockholm subject missing from publication')
    if not matched.paper__sample_type.str.strip().str.upper().eq('TUMOUR').all() or not matched.paper__specimen_type.str.strip().eq('RP').all():
        raise ValueError('Publication does not confirm Stockholm primary RP tumor')
    audit['selected'] = True
    audit['selection_reason'] = 'primary_tumor'
    clinical = cambridge.annotate_clinical(audit.copy().reset_index(drop=True), publication, SERIES)
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title']].copy()
    return clinical, samples, audit
