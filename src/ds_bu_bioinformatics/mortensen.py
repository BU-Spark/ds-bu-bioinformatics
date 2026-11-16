"""Mortensen GSE46602: retain distinct GPL570 probes and source BCR outcomes."""
import numpy as np
import pandas as pd

from .cambridge import read_matrix
from .geo_cohort_pipeline import add_outcomes, decode_bcr

SERIES = 'GSE46602'
PLATFORM = 'GPL570'
STUDY = 'Mortensen / Aarhus'
SNAPSHOT = 'mortensen-v1'
SCALE = 'Deposited RMA-normalized microarray intensity; no added transformation'
MATRIX_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46602/matrix/GSE46602_series_matrix.txt.gz'
PLATFORM_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570&targ=self&form=text&view=full'


def outcome(record):
    return decode_bcr(record.get('geo__bcr'), record.get('geo__bcr_free_time'), {'YES': 1, 'NO': 0})


def clinical_tables(metadata):
    audit = metadata.copy()
    tumor = audit.geo__tissue.eq('prostate tumor') & audit.title.str.fullmatch(r'Prostate cancer cells, PC\d+')
    benign = audit.geo__tissue.eq('benign prostate glands') & audit.title.str.fullmatch(r'Benign prostate glands(?: adjacent to cancer)?, PC\d+')
    if not (tumor | benign).all():
        raise ValueError('Unrecognized Mortensen tissue/title')
    audit['patient_id'] = audit.title.str.rsplit(', ', n=1).str[-1]
    audit['selected'] = tumor
    audit['selection_reason'] = np.where(tumor, 'primary_tumor', 'benign')
    clinical = audit.loc[tumor].reset_index(drop=True)
    if clinical.empty or not clinical.patient_id.is_unique:
        raise ValueError('Empty/duplicate primary Mortensen subjects')
    clinical = add_outcomes(clinical, outcome)
    clinical['bcr_clinical_source'] = 'GEO GSE46602 sample bcr / bcr_free_time'
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title']].copy()
    return clinical, samples, audit


def platform_annotations(path, expression):
    required = ['ID', 'ENTREZ_GENE_ID', 'Gene Symbol', 'RefSeq Transcript ID', 'Species Scientific Name']
    with path.open(encoding='utf-8') as f:
        correct = False
        for line in f:
            correct |= line.rstrip() == '^PLATFORM = GPL570'
            if line.rstrip() == '!platform_table_begin':
                break
        else:
            raise ValueError('GPL570 table absent')
        header = next(f).rstrip('\r\n').split('\t')
        if not correct or not set(required).issubset(header):
            raise ValueError('Unexpected GPL570 schema')
        positions = [header.index(c) for c in required]
        rows = []
        for line in f:
            if line.rstrip() == '!platform_table_end':
                break
            parts = line.rstrip('\r\n').split('\t')
            if len(parts) != len(header):
                raise ValueError('Malformed GPL570 annotation row')
            rows.append([parts[i] for i in positions])
        else:
            raise ValueError('Truncated GPL570 annotation')
    annotations = pd.DataFrame(rows, columns=['feature_id', 'gene_id', 'gene_symbol', 'refseq_original', 'species'])
    if annotations.empty or annotations.feature_id.eq('').any() or not annotations.feature_id.is_unique:
        raise ValueError('Missing/duplicate GPL570 features')
    if not set(expression.index).issubset(annotations.feature_id):
        raise ValueError('Expression features missing from GPL570')
    annotations = annotations.set_index('feature_id').loc[expression.index].reset_index()
    single = annotations.gene_id.str.fullmatch(r'\d+') & annotations.gene_symbol.str.fullmatch(r'[^;/,\s]+') & annotations.species.eq('Homo sapiens')
    pairs = annotations.loc[single, ['gene_id', 'gene_symbol']].drop_duplicates()
    bad_ids = set(pairs.loc[pairs.gene_id.duplicated(False), 'gene_id'])
    bad_symbols = set(pairs.loc[pairs.gene_symbol.duplicated(False), 'gene_symbol'])
    conflict = single & (annotations.gene_id.isin(bad_ids) | annotations.gene_symbol.isin(bad_symbols))
    annotations['mapping_status'] = np.select([conflict, single], ['gene_symbol_identity_conflict', 'mapped'], default='invalid_or_ambiguous_annotation')
    annotations['annotation_eligible'] = annotations.mapping_status.eq('mapped')
    return annotations


def prepare(cache):
    meta, expression, processing = read_matrix(cache.fetch('GSE46602_series_matrix.txt.gz', MATRIX_URL), series=SERIES, platform=PLATFORM)
    clinical, samples, audit = clinical_tables(meta)
    source = expression.loc[:, samples.sample_id]
    annotations = platform_annotations(cache.fetch('GPL570.soft.txt', PLATFORM_URL), source)
    details = {'sample_selection': 'GEO prostate tumor tissue AND Prostate cancer cells, PC* title; benign excluded independently of outcomes',
        'patient_identity': 'PC subject label derived from GEO title; no external registry linkage',
        'annotation': 'Frozen GPL570 single human Entrez/symbol; inconsistent identities quarantined; no fresh MyGene/HGNC validation',
        'outcomes': 'YES/NO BCR and bcr_free_time as deposited; source paper reports postoperative follow-up in months',
        'definition_source': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4629186/',
        'limitations': ['No raw assay QC', 'Deposited platform annotation may be outdated', 'BCR threshold in primary paper PSA >0.2 ng/ml; no independent clinical adjudication']}
    return clinical, samples, audit, source, annotations, {}, processing, details


def summarize(clinical, audit, mapping, extra):
    return {}
