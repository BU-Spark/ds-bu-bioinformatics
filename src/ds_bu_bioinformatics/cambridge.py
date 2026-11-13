"""Strict Cambridge GSE70768 parsing, primary selection and BCR annotation."""
from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .expression import open_text

SERIES = 'GSE70768'
STUDY = 'Cambridge / CamCap'
SNAPSHOT = 'cambridge-v2'
CLINICAL_COHORT = 'cambridge'
SAMPLE_POLICY = 'Select title tumour tissue_robotic radical prostatectomy_* AND sample type Tumour; exclude CRPC and benign independently of outcomes; retain missing clinical outcomes'

MATRIX_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70768/matrix/GSE70768_series_matrix.txt.gz'
PLATFORM_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558&targ=self&form=text&view=full'
MISSING = {'', 'NA', 'N/A', 'UNKNOWN'}


def read_matrix(path: Path, series=SERIES, platform='GPL10558', allow_empty=False):
    headers, characteristics = {}, []
    required = {'!Series_geo_accession', '!Sample_geo_accession', '!Sample_title', '!Sample_platform_id'}
    processing = []
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            if not row:
                continue
            if row == ['!series_matrix_table_begin']:
                break
            if row[0] == '!Sample_characteristics_ch1':
                characteristics.append(row[1:])
            elif row[0] == '!Sample_data_processing':
                processing.append(row[1:])
            elif row[0] in required:
                if row[0] in headers:
                    raise ValueError('Duplicate Illumina GEO header')
                headers[row[0]] = row[1:]
        else:
            raise ValueError('Missing matrix begin')
        if not required.issubset(headers) or headers['!Series_geo_accession'] != [series]:
            raise ValueError(f'Expected {series} metadata')
        gsms = headers['!Sample_geo_accession']
        n = len(gsms)
        values = [headers['!Sample_title'], headers['!Sample_platform_id'], *characteristics, *processing]
        if not n or len(set(gsms)) != n or any(not x for x in gsms) or any(len(v) != n for v in values):
            raise ValueError('Invalid Illumina metadata lengths/accessions')
        if set(headers['!Sample_platform_id']) != {platform} or not processing:
            raise ValueError(f'Expected {platform} and source processing text')
        records = []
        for i, gsm in enumerate(gsms):
            rec = {'sample_id': gsm, 'geo_accession': gsm, 'title': headers['!Sample_title'][i], 'platform_id': platform}
            # Characteristic keys can change between columns within a GEO row.
            for vals in characteristics:
                if not vals[i].strip():
                    continue
                if ':' not in vals[i]:
                    raise ValueError('Malformed Illumina characteristic')
                key, value = vals[i].split(':', 1)
                key = 'geo__' + re.sub('[^a-z0-9]+', '_', key.lower()).strip('_')
                if key in rec:
                    raise ValueError('Duplicate Illumina characteristic')
                rec[key] = value.strip()
            records.append(rec)
        if next(reader, []) != ['ID_REF', *gsms]:
            raise ValueError('Illumina matrix column order differs')
        ids, data = [], []
        for row in reader:
            if row == ['!series_matrix_table_end']:
                break
            if len(row) != n + 1 or not row[0]:
                raise ValueError('Malformed Illumina expression row')
            ids.append(row[0])
            # Only recognized missing tokens are nullable. Invalid numeric text fails.
            parsed = [np.nan if x.strip().upper() in MISSING else float(x) for x in row[1:]]
            if any(not np.isfinite(v) and x.strip().upper() not in MISSING for x, v in zip(row[1:], parsed)):
                raise ValueError('Unrecognized nonfinite Illumina value')
            data.append(parsed)
        else:
            raise ValueError('Truncated Illumina matrix')
    if (not ids and not allow_empty) or len(set(ids)) != len(ids):
        raise ValueError('Empty matrix or duplicate Illumina features')
    expression = pd.DataFrame(data, index=pd.Index(ids, name='feature_id'), columns=gsms, dtype=float)
    if np.isinf(expression.to_numpy()).any():
        raise ValueError('Infinite Illumina expression')
    return pd.DataFrame(records), expression, sorted({x for vals in processing for x in vals})


def clinical_tables(metadata, publication=None):
    audit = metadata.copy()
    if 'geo__sample_type' not in audit:
        raise ValueError('Cambridge sample type absent')
    primary_prefix = 'tumour tissue_robotic radical prostatectomy_'
    crpc_prefix = 'tumour tissue_castrate-resistant prostate cancer (CRPC)_'
    primary = audit.title.str.startswith(primary_prefix) & audit.geo__sample_type.eq('Tumour')
    crpc = audit.title.str.startswith(crpc_prefix) & audit.geo__sample_type.eq('Tumour')
    benign = audit.title.str.startswith(('matched benign tissue_robotic radical prostatectomy_',
                                             'benign tissue_robotic radical prostatectomy_')) & audit.geo__sample_type.eq('Benign')
    if not (primary | crpc | benign).all():
        raise ValueError('Unrecognized Cambridge specimen title/type; review selection')
    audit['selected'] = primary
    audit['selection_reason'] = np.select([primary, crpc, benign], ['primary_tumor', 'crpc', 'benign'], default='unrecognized')
    # This is a GEO-title-derived cohort subject ID, not an independently verified patient registry ID.
    audit['patient_id'] = audit.title.str.rsplit('_', n=1).str[-1]
    clinical = audit.loc[primary].copy().reset_index(drop=True)
    if clinical.empty or clinical.patient_id.eq('').any() or clinical.patient_id.duplicated().any():
        raise ValueError('Missing/duplicate primary Cambridge subject IDs')
    clinical = annotate_clinical(clinical, publication, SERIES)
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title']].copy()
    return clinical, samples, audit


def annotate_clinical(clinical, publication, series):
    if publication is not None:
        clinical = clinical.merge(publication, on='patient_id', how='left', validate='one_to_one', indicator=True, sort=False)
        if not clinical._merge.eq('both').all():
            raise ValueError('Primary GEO subject missing from publication clinical table')
        clinical = clinical.drop(columns='_merge')
        clinical['paper__bcr_status'] = clinical.paper__bcr_status.str.strip()
        clinical['bcr_clinical_source'] = 'PMC4588396 Supplementary Table 2'
    else:
        clinical['bcr_clinical_source'] = 'GEO ' + series
    outcomes = [source_bcr(r) for r in clinical.to_dict('records')]
    clinical['bcr_event'] = pd.array([x[0] for x in outcomes], dtype='Int8')
    clinical['bcr_months'] = pd.array([x[1] for x in outcomes], dtype='Float64')
    clinical['bcr_eligibility_reason'] = [x[2] for x in outcomes]
    clinical['bcr_eligible'] = clinical.bcr_eligibility_reason.eq('eligible')
    clinical['bcr_time_unit'] = 'months'
    return clinical


def source_bcr(record):
    if record.get('bcr_clinical_source') == 'PMC4588396 Supplementary Table 2':
        return decode_bcr(record.get('paper__bcr_status'), record.get('paper__bcr_months'), record.get('paper__followup_months'))
    return decode_bcr(record.get('geo__biochemical_relapse_bcr'), record.get('geo__time_to_bcr_months'), record.get('geo__total_follow_up_months'))


def decode_bcr(status, event_time, followup):
    def missing(x):
        return pd.isna(x) or str(x).strip().upper() in MISSING

    def time(x):
        if missing(x):
            return None
        value = float(x)
        if not np.isfinite(value) or value < 0:
            raise ValueError('Invalid BCR time')
        return value

    event_months, followup_months = time(event_time), time(followup)
    if missing(status):
        return None, None, 'missing_status'
    if status not in {'Y', 'N'}:
        raise ValueError(f'Unknown BCR status: {status!r}')
    event = int(status == 'Y')
    duration = event_months if event else followup_months
    if event and duration is not None and followup_months is not None and duration > followup_months + 1e-8:
        return event, duration, 'event_after_recorded_followup_requires_review'
    reason = 'missing_time' if duration is None else 'zero_time_requires_review' if duration == 0 else 'eligible'
    return event, duration, reason


def map_features(path: Path, expression):
    required = ['ID', 'Entrez_Gene_ID', 'Symbol', 'RefSeq_ID', 'Species']
    records = []
    with path.open(encoding='utf-8') as handle:
        correct = False
        for line in handle:
            correct |= line.rstrip() == '^PLATFORM = GPL10558'
            if line.rstrip() == '!platform_table_begin':
                break
        else:
            raise ValueError('Missing GPL10558 table')
        header = next(handle).rstrip('\r\n').split('\t')
        if not correct or not set(required).issubset(header):
            raise ValueError('Unexpected GPL10558 schema')
        positions = [header.index(k) for k in required]
        for line in handle:
            if line.rstrip() == '!platform_table_end':
                break
            row = line.rstrip('\r\n').split('\t')
            if len(row) != len(header):
                raise ValueError('Malformed GPL10558 row')
            records.append([row[i] for i in positions])
        else:
            raise ValueError('Truncated GPL10558')
    mapping = pd.DataFrame(records, columns=['feature_id', 'gene_id', 'gene_symbol', 'refseq_original', 'species'])
    if mapping.empty or mapping.feature_id.eq('').any() or mapping.feature_id.duplicated().any():
        raise ValueError('Empty/duplicate platform features')
    if not set(expression.index).issubset(mapping.feature_id):
        raise ValueError('Expression features absent from GPL10558')
    mapping = mapping.set_index('feature_id').loc[expression.index].reset_index()
    single = mapping.gene_id.str.fullmatch(r'\d+') & mapping.gene_symbol.str.fullmatch(r'[^;/,\s]+')
    human = mapping.species.eq('Homo sapiens')
    mapping['mapping_status'] = np.where(single & human, 'mapped', 'invalid_or_ambiguous_annotation')
    pairs = mapping.loc[single & human, ['gene_id', 'gene_symbol']].drop_duplicates()
    bad_ids = set(pairs.loc[pairs.gene_id.duplicated(False), 'gene_id'])
    bad_symbols = set(pairs.loc[pairs.gene_symbol.duplicated(False), 'gene_symbol'])
    conflict = single & human & (mapping.gene_id.isin(bad_ids) | mapping.gene_symbol.isin(bad_symbols))
    mapping.loc[conflict, 'mapping_status'] = 'gene_symbol_identity_conflict'
    mapping['missing_selected_values'] = expression.isna().sum(axis=1).to_numpy()
    mapping['expression_status'] = np.where(mapping.missing_selected_values.gt(0), 'missing_selected_expression', 'complete')
    mapping['selected'] = mapping.mapping_status.eq('mapped') & mapping.expression_status.eq('complete')
    mapping['mapping_source'] = 'deposited GPL10558 Entrez_Gene_ID/Symbol; no MyGene refresh'
    genes = mapping.loc[mapping.selected, ['feature_id', 'gene_id', 'gene_symbol', 'refseq_original']].reset_index(drop=True)
    if genes.empty:
        raise ValueError('No finite unambiguously annotated Illumina features')
    return expression.loc[genes.feature_id].copy(), genes, mapping


def validate(expression, clinical, samples, genes, mapping, audit):
    if expression.empty or not np.isfinite(expression.to_numpy()).all():
        raise ValueError('Empty/nonfinite Illumina expression')
    if not expression.index.is_unique or not expression.columns.is_unique:
        raise ValueError('Duplicate Illumina labels')
    if expression.columns.tolist() != samples.sample_id.tolist() or expression.columns.tolist() != clinical.sample_id.tolist():
        raise ValueError('Illumina sample order differs')
    if expression.index.tolist() != genes.feature_id.tolist() or not genes.feature_id.is_unique:
        raise ValueError('Illumina feature order differs')
    if not mapping.feature_id.is_unique or expression.index.tolist() != mapping.loc[mapping.selected, 'feature_id'].tolist():
        raise ValueError('Illumina feature audit differs')
    annotated = mapping.loc[mapping.selected, genes.columns].reset_index(drop=True)
    if not genes.fillna('').astype(str).equals(annotated.fillna('').astype(str)):
        raise ValueError('Illumina gene annotations differ from mapping audit')
    if not samples.equals(clinical[samples.columns]) or samples.patient_id.isna().any() or not samples.patient_id.is_unique:
        raise ValueError('Illumina sample/subject identities differ')
    if not audit.sample_id.is_unique or audit.loc[audit.selected, 'sample_id'].tolist() != samples.sample_id.tolist():
        raise ValueError('Illumina sample audit differs')
    decoded = [source_bcr(r) for r in clinical.to_dict('records')]
    for i, (event, months, reason) in enumerate(decoded):
        row = clinical.iloc[i]
        if (pd.isna(row.bcr_event) != (event is None) or event is not None and row.bcr_event != event
                or pd.isna(row.bcr_months) != (months is None) or months is not None and row.bcr_months != months
                or row.bcr_eligibility_reason != reason or bool(row.bcr_eligible) != (reason == 'eligible')):
            raise ValueError('Illumina BCR fields disagree with source')
    eligible = clinical.bcr_eligible
    return {'structurally_valid': True, 'counts': {'primary_profiles': len(samples), 'retained_feature_rows': len(genes),
        'unique_genes': int(genes.gene_id.nunique()), 'bcr_events': int(clinical.bcr_event.eq(1).sum()),
        'non_events': int(clinical.bcr_event.eq(0).sum()), 'eligible_bcr_profiles': int(eligible.sum()),
        'eligible_bcr_events': int(clinical.loc[eligible, 'bcr_event'].eq(1).sum()),
        'complete_positive_time_bcr_profiles': int((clinical.bcr_event.notna() & clinical.bcr_months.gt(0)).sum())},
        'geo_profiles': len(audit), 'source_features': len(mapping),
        'sample_dispositions': audit.selection_reason.value_counts().to_dict(),
        'mapping_dispositions': mapping.mapping_status.value_counts().to_dict(),
        'features_with_missing_selected_expression': int(mapping.missing_selected_values.gt(0).sum()),
        'missing_selected_expression_cells': int(mapping.missing_selected_values.sum()),
        'mapped_features_excluded_for_missing_expression': int((mapping.mapping_status.eq('mapped') & mapping.missing_selected_values.gt(0)).sum()),
        'bcr_eligibility_reasons': clinical.bcr_eligibility_reason.value_counts().to_dict(),
        'clinical_missingness': clinical.isna().sum().to_dict(),
        'expression_range': [float(expression.to_numpy().min()), float(expression.to_numpy().max())]}
