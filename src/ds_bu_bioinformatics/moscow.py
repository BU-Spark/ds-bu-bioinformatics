"""Moscow/EMB-RAS GSE229904 mRNA counts, assay selection and deposited BCR labels."""
import re

import numpy as np
import pandas as pd

from .gene_identifiers import annotate
from .geo_cohort_pipeline import add_outcomes, decode_bcr
from .geo_io import align_columns, read_geo, read_tsv

SERIES = 'GSE229904'
PLATFORM = 'GPL18573'
STUDY = 'Moscow / EMB-RAS'
SNAPSHOT = 'moscow-v1'
SCALE = 'Unmodified deposited mRNA integer counts; raw exports; explicit log-CPM only during module scoring'
BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229904'
MATRIX_URL = BASE + '/matrix/GSE229904_series_matrix.txt.gz'
COUNTS_URL = BASE + '/suppl/GSE229904_mRNA_counts.tsv.gz'


def outcome(record):
    return decode_bcr(record.get('geo__biochemical_recurrence'), None, {'bcr': 1, 'brf': 0})


ENDPOINTS = {'bcr': outcome}


def identifier_scope(ids):
    if all(re.fullmatch(r'ENSG\d+(?:\.\d+)?', x) for x in ids):
        return 'ensembl.gene'
    if all(re.fullmatch(r'\d+(?:_at)?', x) for x in ids):
        return 'entrezgene'
    if all(re.fullmatch(r'[A-Za-z][A-Za-z0-9_.-]*', x) for x in ids):
        return 'symbol'
    raise ValueError('Mixed/unsupported source gene identifiers; inspect source annotation before scoring')


def clinical_tables(metadata, present):
    audit = metadata.copy()
    if audit.empty or not audit.sample_id.is_unique or not audit.platform_id.eq(PLATFORM).all() or 'geo__tissue' not in audit:
        raise ValueError('Unexpected Moscow assay metadata')
    rna = audit.title.str.fullmatch(r'\d+[TN] RNA')
    mirna = audit.title.str.fullmatch(r'\d+[TN] miRNA')
    if not (rna | mirna).all():
        raise ValueError('Unrecognized Moscow assay type/title')
    primary = audit.geo__tissue.eq('primary tumor') & audit.title.str.match(r'^\d+T ')
    benign = audit.geo__tissue.eq('normal prostate tissue') & audit.title.str.match(r'^\d+N ')
    if not (primary | benign).all():
        raise ValueError('Moscow specimen title and tissue disagree')
    audit['patient_id'] = SERIES + ':' + audit.title.str.extract(r'^(\d+)', expand=False)
    has_counts = audit.sample_id.isin(present)
    if (has_counts & ~rna).any():
        raise ValueError('miRNA assay was linked to an mRNA count column')
    audit['selected'] = rna & primary & has_counts
    audit['selection_reason'] = np.select([~rna, benign, rna & primary & ~has_counts, audit.selected],
        ['miRNA_assay', 'normal_tissue', 'mRNA_primary_absent_from_counts', 'primary_mRNA'], default='unrecognized')
    clinical = add_outcomes(audit.loc[audit.selected].reset_index(drop=True), outcome)
    if clinical.empty or not clinical.patient_id.is_unique:
        raise ValueError('Missing/duplicate Moscow primary mRNA patients')
    clinical['bcr_clinical_source'] = 'GEO biochemical recurrence bcr/brf class; individual durations not verified'
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title']].copy()
    return clinical, samples, audit


def prepare(cache):
    meta, _, processing = read_geo(cache.fetch(SERIES + '_series_matrix.txt.gz', MATRIX_URL),
        series=SERIES, platforms={PLATFORM}, metadata_only=True)
    raw = read_tsv(cache.fetch('GSE229904_mRNA_counts.tsv.gz', COUNTS_URL))
    rna = meta.loc[meta.title.str.fullmatch(r'\d+[TN] RNA')]
    raw, linkage = align_columns(raw, rna,
        lambda r: [r['sample_id'], r['title'], r['title'].removesuffix(' RNA')])
    clinical, samples, audit = clinical_tables(meta, set(raw.columns))
    source = raw.loc[:, samples.sample_id].copy()
    values = source.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError('Moscow raw counts must be complete nonnegative integers')
    totals = source.sum(axis=0)
    if not np.isfinite(totals).all() or not totals.gt(0).all():
        raise ValueError('Zero/nonfinite mRNA library size')
    scope = identifier_scope(source.index)
    genes, candidates = annotate(cache, source.index.tolist(), scope=scope)
    details = {'sample_selection': 'Primary tumor mRNA only; exclude normal/miRNA assays; exact GSM/title/core specimen key joins; all deposited assays audited',
        'patient_identity': 'Numeric source specimen title label, before T/N suffix; repeated selected patient fails',
        'annotation': 'Exact source gene namespace mapped to human Entrez/symbol with cached responses; original gene rows retained',
        'expression_kind': 'raw_counts', 'library_sizes': {str(k): float(v) for k, v in totals.items()},
        'endpoints': {'bcr': {'definition': 'Deposited bcr/brf biochemical recurrence class', 'duration_available': False, 'unit': 'months'}},
        'study_source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229904',
        'limitations': ['Full-source schema/count verification pending; synthetic tests only',
            'Published patient totals are not assumed to equal mRNA profiles with usable BCR metadata',
            'No verified individual BCR/censoring times; survival requires a linked supplement',
            'Raw counts remain unchanged; log-CPM scoring is basic library-size normalization and does not resolve batch/composition effects']}
    return clinical, samples, audit, source, genes, {'expression_linkage.csv': linkage,
        'gene_mapping_candidates.csv.gz': candidates}, processing, details


def summarize(clinical, audit, mapping, extra):
    return {'survival_ready': bool(clinical.bcr_eligible.any()),
        'primary_mRNA_profiles_absent_from_counts': int(audit.selection_reason.eq('mRNA_primary_absent_from_counts').sum())}
