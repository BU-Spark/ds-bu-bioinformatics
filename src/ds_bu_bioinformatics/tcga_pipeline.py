"""Download, audit, and export TCGA-PRAD (cBioPortal prad_tcga_pan_can_atlas_2018) RNA-seq + clinical data.

Output layout and column conventions follow the Taylor pipeline (see docs/data_dictionary.md):
expression.csv.gz (feature_id x sample_id), genes.csv, samples.csv, clinical.csv (sample__/patient__
source fields + derived endpoint columns), sample_audit.csv, feature_audit.csv, metadata.json,
validation.json, checksums.sha256, and an AnnData file in a separate folder.

TCGA-specific decisions (CampLab client meetings, Oct-Nov 2026; team meeting 11/06):
  * Expression is log2(RSEM + 1) after clipping small negative batch-corrected values to 0. Raw RSEM
    values are kept in the AnnData layer 'rsem'. Every source row with complete values is kept;
    rows sharing an Entrez ID or symbol stay separate (no collapsing).
  * Primary endpoint: progression-free interval (PFI), cBioPortal PFS_STATUS / PFS_MONTHS.
    Secondary: OS, DSS, DFI. Each uses the same nullable event/time/eligibility pattern as Taylor DFS.
  * Sample selection: barcode type 01 (primary tumor), one sample per patient (lowest barcode),
    then exclusion of expression QC outliers (|z| > threshold on median expression or genes detected).
  * Gleason comes from the TCGA Firehose Legacy study (prad_tcga), matched on patient ID.
"""
from __future__ import annotations

import argparse
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .anndata_export import write_anndata
from .clinical import paged, pivot_clinical
from .expression import write_csv
from .sources import SourceCache, sha256, write_json
from .taylor_pipeline import check_destination, publish

STUDY = 'prad_tcga_pan_can_atlas_2018'
LEGACY_STUDY = 'prad_tcga'
CBIO_API = 'https://www.cbioportal.org/api'
EXPR_NAME = 'data_mrna_seq_v2_rsem.txt'
EXPR_URL = f'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/{STUDY}/{EXPR_NAME}'
PLATFORM_ID = 'IlluminaHiSeq_RNASeqV2'
DETECTED_LOG2 = 1.0   # a feature counts as detected in a sample above log2(RSEM + 1) = 1
DECIMALS = 4          # stored precision of log2 expression
PAGE_SIZE = 10000
BARCODE = r'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}[A-Z]?'
MISSING_TOKENS = {'', 'NA', 'N/A', '[Not Available]', '[Not Applicable]', '[Unknown]', '[Discrepancy]'}

# Endpoint prefix -> (status attribute, months attribute, explicit status code map). Unknown codes fail.
ENDPOINTS = {
    'pfi': ('PFS_STATUS', 'PFS_MONTHS', {'0:CENSORED': 0, '1:PROGRESSION': 1}),
    'os': ('OS_STATUS', 'OS_MONTHS', {'0:LIVING': 0, '1:DECEASED': 1}),
    'dss': ('DSS_STATUS', 'DSS_MONTHS', {'0:ALIVE OR DEAD TUMOR FREE': 0, '1:DEAD WITH TUMOR': 1}),
    'dfi': ('DFS_STATUS', 'DFS_MONTHS', {'0:DiseaseFree': 0, '1:Recurred/Progressed': 1}),
}
PRIMARY_ENDPOINT = 'pfi'


# ----------------------------------------------------------------------------- expression
def read_rsem(path: Path):
    """Parse the cBioPortal RSEM matrix. Returns (raw values features x samples, feature table)."""
    with open(path, encoding='utf-8') as handle:
        header = handle.readline().rstrip('\r\n').split('\t')
    if header[:2] != ['Hugo_Symbol', 'Entrez_Gene_Id']:
        raise ValueError('Expected Hugo_Symbol and Entrez_Gene_Id as the first RSEM columns')
    samples = pd.Index(header[2:])
    # Check the raw header: pandas would silently rename a repeated column to <name>.1
    if samples.empty or not samples.is_unique or not samples.str.fullmatch(BARCODE).all():
        raise ValueError('RSEM sample columns must be unique TCGA barcodes')
    frame = pd.read_csv(path, sep='\t', dtype={'Hugo_Symbol': 'string'}, na_values=['NA', 'NaN', ''],
                        keep_default_na=False, low_memory=False)
    if frame.columns.tolist() != header:
        raise ValueError('RSEM header changed while parsing')
    entrez = pd.to_numeric(frame['Entrez_Gene_Id'], errors='coerce')
    if entrez.isna().any() or (entrez <= 0).any() or (entrez % 1 != 0).any():
        raise ValueError('Every RSEM row needs a positive integer Entrez_Gene_Id')
    try:
        values = frame[samples].astype(float)
    except ValueError as error:
        raise ValueError(f'Non-numeric RSEM value: {error}') from error
    gene_id = entrez.astype('int64').astype(str).to_numpy(dtype=object)
    source_row = np.arange(1, len(frame) + 1)
    repeated = pd.Series(gene_id).duplicated(keep=False).to_numpy()
    # Entrez ID is the row key; the few repeated IDs get the 1-based source data row appended.
    feature_id = np.array([f'{g}_row{r}' if rep else g for g, r, rep in zip(gene_id, source_row, repeated)], dtype=object)
    symbol = frame['Hugo_Symbol'].str.strip().reset_index(drop=True)
    features = pd.DataFrame({'feature_id': feature_id, 'source_row': source_row, 'gene_id': gene_id,
                             'gene_symbol': symbol.mask(symbol.isin(MISSING_TOKENS))})
    features = features.astype({'feature_id': 'string', 'gene_id': 'string', 'gene_symbol': 'string'})
    values.index = pd.Index(features.feature_id.astype(str), name='feature_id')
    values.columns = pd.Index(samples.astype(str), name=None)
    if not values.index.is_unique:
        raise ValueError('Duplicate feature IDs after assignment')
    return values, features


def clean_expression(raw: pd.DataFrame, features: pd.DataFrame):
    """Drop rows with missing values, clip negatives to 0, log2(x + 1). Returns (expr, raw kept, audit, stats)."""
    if np.isinf(raw.to_numpy()).any():
        raise ValueError('Infinite RSEM values')
    n_missing = raw.isna().sum(axis=1)
    audit = features.copy()
    audit['missing_samples'] = n_missing.to_numpy()
    audit['selected'] = audit.missing_samples.eq(0)
    audit['reason'] = np.where(audit.selected, 'selected', 'missing_values')
    kept = raw.loc[audit.loc[audit.selected, 'feature_id'].astype(str)]
    if kept.empty:
        raise ValueError('No features without missing values')
    negative = int((kept < 0).sum().sum())
    expression = np.log2(kept.clip(lower=0) + 1).round(DECIMALS)
    stats = {'source_features': len(raw), 'source_profiles': raw.shape[1],
             'features_missing_values': int((~audit.selected).sum()),
             'samples_with_missing_values': int(raw.isna().any(axis=0).sum()),
             'negative_rsem_values_clipped': negative}
    return expression, kept, audit, stats


# ----------------------------------------------------------------------------- samples
def select_samples(expression: pd.DataFrame, z_threshold: float = 3.0) -> pd.DataFrame:
    """Sample audit for every expression column: barcode rule, one per patient, QC outliers."""
    audit = pd.DataFrame({'sample_id': expression.columns.astype(str)})
    audit['patient_id'] = audit.sample_id.str[:12]
    audit['sample_type_code'] = audit.sample_id.str[13:15]
    audit['selected'] = audit.sample_type_code.eq('01')
    audit['selection_reason'] = np.where(audit.selected, 'primary_tumor', 'not_primary_tumor')
    extra = audit[audit.selected].sort_values('sample_id').duplicated('patient_id')
    extra = extra[extra].index
    audit.loc[extra, ['selected', 'selection_reason']] = [False, 'extra_primary_for_patient']

    kept = audit.loc[audit.selected, 'sample_id']
    median = expression[kept].median()
    detected = (expression[kept] > DETECTED_LOG2).sum()

    def z(series):
        sd = series.std(ddof=0)
        return (series - series.mean()) / sd if sd > 0 else series * 0.0

    outlier = (z(median).abs() > z_threshold) | (z(detected).abs() > z_threshold)
    audit['median_log2'] = audit.sample_id.map(expression.median().round(DECIMALS))
    audit['genes_detected'] = audit.sample_id.map((expression > DETECTED_LOG2).sum()).astype('int64')
    is_outlier = audit.sample_id.isin(outlier[outlier].index)
    audit.loc[is_outlier, ['selected', 'selection_reason']] = [False, 'qc_outlier']
    return audit


# ----------------------------------------------------------------------------- clinical
def _missing(value) -> bool:
    return pd.isna(value) or str(value).strip() in MISSING_TOKENS


def decode_endpoint(status, raw_time, codes):
    """Return (event, months, eligibility_reason). Missing stays missing; unknown/invalid values fail."""
    missing_status, missing_time = _missing(status), _missing(raw_time)
    if not missing_status and status not in codes:
        raise ValueError(f'Unknown endpoint status: {status!r}')
    event = None if missing_status else codes[status]
    months = None if missing_time else float(raw_time)
    if months is not None and (not np.isfinite(months) or months < 0):
        raise ValueError(f'Invalid endpoint time: {raw_time!r}')
    reason = ('missing_status' if missing_status else 'missing_time' if missing_time
              else 'zero_time_requires_review' if months == 0 else 'eligible')
    return event, months, reason


def _column(frame: pd.DataFrame, name: str) -> pd.Series:
    return frame[name] if name in frame else pd.Series(None, index=frame.index, dtype=object)


def add_endpoints(clinical: pd.DataFrame) -> pd.DataFrame:
    for prefix, (status_col, time_col, codes) in ENDPOINTS.items():
        status = _column(clinical, f'patient__{status_col}')
        times = _column(clinical, f'patient__{time_col}')
        decoded = [decode_endpoint(s, t, codes) for s, t in zip(status, times, strict=True)]
        clinical[f'{prefix}_event'] = pd.array([d[0] for d in decoded], dtype='Int8')
        clinical[f'{prefix}_months'] = pd.array([d[1] for d in decoded], dtype='Float64')
        clinical[f'{prefix}_eligibility_reason'] = [d[2] for d in decoded]
        clinical[f'{prefix}_eligible'] = clinical[f'{prefix}_eligibility_reason'].eq('eligible')
        clinical[f'{prefix}_status'] = status.map(lambda x: None if _missing(x) else str(x).split(':', 1)[-1])
        clinical[f'{prefix}_time_unit'] = 'months'
    return clinical


def match_clinical(audit, patient_records, sample_records, inventory_records, legacy_records):
    """Join selected samples to cBioPortal records. Returns (clinical, samples, updated audit, stats)."""
    patient, p_dups = pivot_clinical(patient_records, 'PATIENT', study=STUDY)
    sample, s_dups = pivot_clinical(sample_records, 'SAMPLE', study=STUDY)
    gleason = [r for r in legacy_records if 'GLEASON' in str(r.get('clinicalAttributeId', '')).upper()]
    legacy, l_dups = pivot_clinical(gleason, 'PATIENT', study=LEGACY_STUDY) if gleason else (
        pd.DataFrame({'patientId': pd.Series(dtype=str)}), 0)
    inventory = pd.DataFrame(inventory_records)
    if not {'sampleId', 'patientId'}.issubset(inventory.columns):
        raise ValueError('Missing cBio sample inventory identifiers')
    inventory = inventory[['sampleId', 'patientId', *(['sampleType'] if 'sampleType' in inventory else [])]]
    if inventory[['sampleId', 'patientId']].isna().any().any() or inventory.sampleId.duplicated().any():
        raise ValueError('Duplicate/missing cBio sample inventory IDs')
    relation = sample[['sampleId', 'patientId']].merge(inventory, on='sampleId', how='left',
                                                       suffixes=('', '_inventory'), validate='one_to_one')
    if not relation.patientId.eq(relation.patientId_inventory).all():
        raise ValueError('Clinical patient mapping disagrees with sample inventory')

    audit = audit.merge(inventory.rename(columns={'sampleId': 'cbio_sample_id', 'patientId': 'cbio_patient_id',
                                                  'sampleType': 'cbio_sample_type'}),
                        left_on='sample_id', right_on='cbio_sample_id', how='left', validate='one_to_one')
    mismatch = audit.cbio_patient_id.notna() & audit.cbio_patient_id.ne(audit.patient_id)
    if mismatch.any():
        raise ValueError('Barcode patient ID disagrees with cBio inventory: ' + ','.join(audit.loc[mismatch, 'sample_id']))
    audit = audit.drop(columns='cbio_patient_id')
    no_inventory = audit.selected & audit.cbio_sample_id.isna()
    audit.loc[no_inventory, ['selected', 'selection_reason']] = [False, 'not_in_cbio_inventory']
    no_clinical = audit.selected & (~audit.sample_id.isin(sample.sampleId) | ~audit.patient_id.isin(patient.patientId))
    audit.loc[no_clinical, ['selected', 'selection_reason']] = [False, 'no_clinical_record']

    selected = audit.loc[audit.selected, ['sample_id', 'cbio_sample_id', 'patient_id', 'sample_type_code']]
    if selected.empty:
        raise ValueError('No selected primary tumors with clinical records')
    if selected.patient_id.duplicated().any():
        raise ValueError('Multiple selected samples per patient')
    sample = sample.rename(columns={c: 'sample__' + c for c in sample.columns if c not in {'sampleId', 'patientId'}})
    patient = patient.rename(columns={c: 'patient__' + c for c in patient.columns if c != 'patientId'})
    legacy = legacy.rename(columns={c: 'legacy_patient__' + c for c in legacy.columns if c != 'patientId'})
    clinical = (selected.merge(sample.drop(columns='patientId'), left_on='sample_id', right_on='sampleId',
                               how='left', validate='one_to_one').drop(columns='sampleId')
                .merge(patient, left_on='patient_id', right_on='patientId', how='left', validate='many_to_one')
                .drop(columns='patientId')
                .merge(legacy, left_on='patient_id', right_on='patientId', how='left', validate='many_to_one')
                .drop(columns='patientId'))
    if 'sample__SAMPLE_TYPE' in clinical and not clinical['sample__SAMPLE_TYPE'].dropna().eq('Primary').all():
        raise ValueError('Barcode and cBio primary tumor labels disagree')
    missing_endpoint = [c for c in ENDPOINTS[PRIMARY_ENDPOINT][:2] if f'patient__{c}' not in clinical]
    if missing_endpoint:
        raise ValueError(f'Primary endpoint fields absent: {missing_endpoint}')
    clinical['age'] = pd.to_numeric(_column(clinical, 'patient__AGE'), errors='coerce').astype('Float64')
    t_stage = _column(clinical, 'patient__PATH_T_STAGE').astype('string')
    clinical['path_t_group'] = t_stage.str.extract(r'(T\d)', expand=False)
    score = pd.to_numeric(_column(clinical, 'legacy_patient__GLEASON_SCORE'), errors='coerce')
    clinical['gleason_score'] = score.astype('Float64').astype('Int8')
    clinical = add_endpoints(clinical)
    samples = clinical[['sample_id', 'cbio_sample_id', 'patient_id', 'sample_type_code']].copy()
    samples['platform_id'] = PLATFORM_ID
    if 'cbio_sample_type' in audit:
        samples['cbio_sample_type'] = samples.sample_id.map(audit.set_index('sample_id').cbio_sample_type)
    stats = {'identical_patient_records_collapsed': p_dups, 'identical_sample_records_collapsed': s_dups,
             'identical_legacy_records_collapsed': l_dups}
    return clinical, samples, audit, stats


# ----------------------------------------------------------------------------- validation
def validate(expression, clinical, samples, genes, feature_audit, audit):
    if expression.empty or not np.isfinite(expression.to_numpy()).all():
        raise ValueError('Empty or non-finite final expression')
    if not expression.index.is_unique or not expression.columns.is_unique:
        raise ValueError('Duplicate expression labels')
    if expression.columns.tolist() != samples.sample_id.tolist() or expression.columns.tolist() != clinical.sample_id.tolist():
        raise ValueError('Expression/sample/clinical order differs')
    if expression.index.astype(str).tolist() != genes.feature_id.astype(str).tolist():
        raise ValueError('Expression/gene annotation order differs')
    if not genes.feature_id.is_unique or set(expression.index.astype(str)) != set(
            feature_audit.loc[feature_audit.selected, 'feature_id'].astype(str)):
        raise ValueError('Selected features lost or duplicated')
    if samples.patient_id.isna().any() or not samples.patient_id.is_unique:
        raise ValueError('Missing/duplicate selected patients')
    identity = ['sample_id', 'cbio_sample_id', 'patient_id']
    if not samples[identity].astype(str).equals(clinical[identity].astype(str)):
        raise ValueError('Clinical/sample identity mismatch')
    if audit.sample_id.duplicated().any() or set(audit.loc[audit.selected, 'sample_id']) != set(samples.sample_id):
        raise ValueError('Sample audit selection mismatch')
    endpoint_counts = {}
    for prefix in ENDPOINTS:
        event, months, eligible = clinical[f'{prefix}_event'], clinical[f'{prefix}_months'], clinical[f'{prefix}_eligible']
        if not event.dropna().isin([0, 1]).all():
            raise ValueError(f'Invalid {prefix} event codes')
        if clinical.loc[eligible, [f'{prefix}_event', f'{prefix}_months']].isna().any().any() or not months[eligible].gt(0).all():
            raise ValueError(f'Invalid eligible {prefix} times/status')
        endpoint_counts[prefix] = {'eligible': int(eligible.sum()), 'events': int(event[eligible].eq(1).sum()),
                                   'censored': int(event[eligible].eq(0).sum()),
                                   'median_observed_months': float(months[eligible].median()) if eligible.any() else None}
    p = PRIMARY_ENDPOINT
    counts = {'primary_profiles': len(samples), f'{p}_events': int(clinical[f'{p}_event'].eq(1).sum()),
              f'{p}_censored': int(clinical[f'{p}_event'].eq(0).sum()), 'gene_rows': int(genes.gene_id.nunique())}
    return {'structurally_valid': True, 'counts': counts, 'primary_endpoint': p, 'endpoints': endpoint_counts,
            'retained_feature_rows': len(genes), 'source_profiles': len(audit),
            'source_features': int(feature_audit.feature_id.nunique()),
            'feature_dispositions': feature_audit.reason.value_counts().to_dict(),
            'sample_dispositions': audit.selection_reason.value_counts().to_dict(),
            f'eligible_{p}_profiles': endpoint_counts[p]['eligible'],
            'clinical_missingness': {k: int(v) for k, v in clinical.isna().sum().items()},
            'expression_range': [float(expression.to_numpy().min()), float(expression.to_numpy().max())]}


def read_export(path: Path):
    expression = pd.read_csv(path / 'expression.csv.gz', index_col=0, dtype={'feature_id': 'string'},
                             float_precision='round_trip')
    expression.index = expression.index.astype(str)
    ids = {c: 'string' for c in ['sample_id', 'cbio_sample_id', 'patient_id', 'sample_type_code']}
    clinical = pd.read_csv(path / 'clinical.csv', dtype=ids, low_memory=False)
    samples = pd.read_csv(path / 'samples.csv', dtype='string')
    genes = pd.read_csv(path / 'genes.csv', dtype={'feature_id': 'string', 'gene_id': 'string', 'gene_symbol': 'string'})
    feature_audit = pd.read_csv(path / 'feature_audit.csv', dtype={'feature_id': 'string', 'gene_id': 'string'})
    audit = pd.read_csv(path / 'sample_audit.csv', dtype=ids)
    return expression, clinical, samples, genes, feature_audit, audit


def validate_export(path: Path, expected_expression=None):
    expression, clinical, samples, genes, feature_audit, audit = read_export(path)
    if expected_expression is not None:
        if (expression.index.tolist() != expected_expression.index.astype(str).tolist()
                or expression.columns.tolist() != expected_expression.columns.tolist()):
            raise ValueError('Export changed expression labels')
        if not np.array_equal(expression.to_numpy(), expected_expression.to_numpy()):
            raise ValueError('Export changed expression values')
    return validate(expression, clinical, samples, genes, feature_audit, audit)


# ----------------------------------------------------------------------------- module coverage
def module_coverage(module_file: Path, genes: pd.DataFrame, feature_audit: pd.DataFrame) -> pd.DataFrame:
    """Share of each module's genes measured. File: CSV/TSV with a gene/symbol column and a module column."""
    mods = pd.read_csv(module_file, sep=None, engine='python', dtype=str)
    gene_col = next((c for c in mods.columns if c.lower() in ('gene', 'genes', 'symbol', 'gene_symbol', 'hugo_symbol')), None)
    if gene_col is None or len(mods.columns) < 2:
        raise ValueError('Module file needs a gene/symbol column and a module column')
    module_col = next(c for c in mods.columns if c != gene_col)
    symbol = mods[gene_col].str.strip()
    mods['found'] = symbol.isin(set(genes.gene_symbol.dropna()))
    mods['lost_to_missing_values'] = symbol.isin(set(feature_audit.loc[~feature_audit.selected, 'gene_symbol'].dropna())) & ~mods.found
    coverage = mods.groupby(module_col).agg(genes=('found', 'size'), found=('found', 'sum'),
                                            lost_to_missing_values=('lost_to_missing_values', 'sum'))
    coverage['pct_found'] = (100 * coverage.found / coverage.genes).round(1)
    return coverage.reset_index()


# ----------------------------------------------------------------------------- reporting
def qc_summary(checks: dict, stats: dict, z_threshold: float) -> str:
    p = PRIMARY_ENDPOINT
    ep = checks['endpoints'][p]
    lines = ['# TCGA-PRAD QC summary', '',
             '## Expression',
             f"- Source matrix: {stats['source_features']:,} features x {stats['source_profiles']} samples",
             f"- Removed {stats['features_missing_values']:,} features with missing values "
             f"(in {stats['samples_with_missing_values']} samples; PanCanAtlas batch-correction gaps); all samples kept",
             f"- Clipped {stats['negative_rsem_values_clipped']:,} negative batch-corrected RSEM values to 0 before log2(x + 1)",
             f"- Final: {checks['retained_feature_rows']:,} feature rows, {checks['counts']['gene_rows']:,} unique Entrez genes; "
             'rows sharing a gene ID or symbol are kept separately',
             '', '## Samples']
    lines += [f'- {reason}: {n}' for reason, n in sorted(checks['sample_dispositions'].items())]
    lines += [f'- QC outlier rule: |z| > {z_threshold:g} on median log2 expression or genes detected',
              '', '## Endpoints (eligible = status and positive time present)']
    for prefix, values in checks['endpoints'].items():
        median = values['median_observed_months']
        lines.append(f"- {prefix.upper()}: {values['eligible']} eligible, {values['events']} events, "
                     f"{values['censored']} censored; median observed time "
                     f"{'n/a' if median is None else f'{median:.1f}'} months")
    lines += ['', '## Output',
              f"- {checks['counts']['primary_profiles']} samples (1 per patient, primary tumor) x "
              f"{checks['retained_feature_rows']:,} features; {ep['events']} {p.upper()} events",
              f"- Export round trip valid: {checks.get('export_round_trip_valid')}; "
              f"AnnData round trip valid: {checks.get('anndata_round_trip_valid')}"]
    missing = {k: v for k, v in checks['clinical_missingness'].items()
               if v and not k.startswith(('patient__', 'sample__', 'legacy_patient__'))}
    if missing:
        lines += ['', '## Missing derived clinical values (samples)'] + [f'- {k}: {v}' for k, v in missing.items()]
    return '\n'.join(lines) + '\n'


def build_metadata(cache, attributes, z_threshold, stats, module_file, clinical, samples, genes, feature_audit, audit):
    attribute_map = {a['clinicalAttributeId']: a for a in attributes if a.get('patientAttribute')}
    needed = [c for status, time, _ in ENDPOINTS.values() for c in (status, time)]
    metadata = {
        'schema_version': 1, 'study': 'TCGA-PRAD', 'cbio_study': STUDY, 'legacy_cbio_study': LEGACY_STUDY,
        'snapshot': cache.root.name, 'source_manifest_sha256': sha256(cache.manifest_path),
        'expression': {
            'orientation': 'features_by_samples', 'source_file': EXPR_NAME, 'source_url': EXPR_URL,
            'feature_key': 'Entrez gene ID; when an ID repeats in the source, <entrez>_row<source data row>',
            'source_scale': 'RSEM normalized counts, PanCanAtlas batch-corrected (EB++); can be slightly negative',
            'transformations': f'clip negative values to 0, then log2(x + 1), rounded to {DECIMALS} decimals',
            'raw_values': "AnnData layer 'rsem': deposited RSEM values for retained features, unclipped",
            'aggregation': 'none; rows sharing an Entrez ID or symbol are kept separately',
            'feature_filter': 'features with a missing value in any sample are removed (all samples kept)',
            'negative_rsem_values_clipped': stats['negative_rsem_values_clipped']},
        'endpoints': {prefix: {'status': f'patient__{s}', 'months': f'patient__{t}', 'status_codes': codes,
                               'definitions': {k: attribute_map[k] for k in (s, t) if k in attribute_map}}
                      for prefix, (s, t, codes) in ENDPOINTS.items()},
        'primary_endpoint': PRIMARY_ENDPOINT,
        'endpoint_policy': {'unit': 'months', 'zero_time_policy': 'retain but exclude from eligibility pending review',
                            'missing_policy': 'retain nullable event/time and eligibility reason; never assume censored',
                            'source': 'TCGA Pan-Cancer Clinical Data Resource endpoints as exposed by cBioPortal; '
                                      'PFS_* is the CDR progression-free interval (PFI), DFS_* the disease-free interval (DFI)'},
        'sample_policy': {'barcode': 'characters 1-12 = patient, 14-15 = sample type; keep 01 (primary tumor)',
                          'one_per_patient': 'lowest sample barcode kept (outcome-independent)',
                          'qc_outliers': f'exclude |z| > {z_threshold:g} on median log2 expression or genes detected '
                                         f'(log2(RSEM + 1) > {DETECTED_LOG2:g})',
                          'identity': 'barcode patient must equal cBio inventory patient; cBio SAMPLE_TYPE must be Primary'},
        'clinical_field_provenance': {
            'patient__*': f'{STUDY} PATIENT clinical records, pivoted by patientId',
            'sample__*': f'{STUDY} SAMPLE clinical records, pivoted by sampleId',
            'legacy_patient__GLEASON_*': f'{LEGACY_STUDY} (Firehose Legacy) PATIENT Gleason records, matched on patientId',
            'age': 'numeric patient__AGE (years at diagnosis)',
            'path_t_group': 'T category (T2/T3/T4) extracted from patient__PATH_T_STAGE',
            'gleason_score': 'numeric legacy_patient__GLEASON_SCORE',
            '<endpoint>_*': 'derived from the endpoint status/months fields listed under endpoints'},
        'clinical_attributes': attributes,
        'limitations': [
            'PanCanAtlas batch correction leaves small negative RSEM values; they are clipped to 0 before log2',
            'Features with batch-correction gaps (missing values) are removed rather than imputed',
            'PFI is not relabeled BCR or DFS; DFI is incomplete in PRAD',
            'Gleason comes from a different cBioPortal study (Firehose Legacy) matched on patient ID'],
        'code_sha256': {p.name: sha256(p) for p in sorted(Path(__file__).parent.glob('*.py'))}}
    if module_file is not None:
        metadata['module_file'] = {'name': Path(module_file).name, 'sha256': sha256(Path(module_file))}
    missing_defs = [k for k in needed if k not in attribute_map]
    if missing_defs:
        metadata['limitations'].append(f'No cBioPortal attribute definition for: {missing_defs}')
    metadata['output_schemas'] = {name: {c: {'dtype': str(f[c].dtype), 'nullable': bool(f[c].isna().any())} for c in f.columns}
                                  for name, f in [('clinical.csv', clinical), ('samples.csv', samples), ('genes.csv', genes),
                                                  ('feature_audit.csv', feature_audit), ('sample_audit.csv', audit)]}
    return metadata


# ----------------------------------------------------------------------------- build
def acquire(cache):
    cache.verify_all()
    base = f'{CBIO_API}/studies/{STUDY}'
    legacy_base = f'{CBIO_API}/studies/{LEGACY_STUDY}'
    expression_path = cache.fetch(EXPR_NAME, EXPR_URL)
    attributes = cache.json('cbio_clinical_attributes.json', f'{base}/clinical-attributes')
    summary = {'projection': 'SUMMARY'}
    patient = paged(cache, 'clinical-data', 'cbio_patient_clinical', {'clinicalDataType': 'PATIENT', **summary},
                    page_size=PAGE_SIZE, base=base)
    sample = paged(cache, 'clinical-data', 'cbio_sample_clinical', {'clinicalDataType': 'SAMPLE', **summary},
                   page_size=PAGE_SIZE, base=base)
    inventory = paged(cache, 'samples', 'cbio_samples', summary, page_size=PAGE_SIZE, base=base)
    legacy = paged(cache, 'clinical-data', 'legacy_patient_clinical', {'clinicalDataType': 'PATIENT', **summary},
                   page_size=PAGE_SIZE, base=legacy_base)
    return expression_path, attributes, patient, sample, inventory, legacy


def build(cache, output_root: Path, anndata_root: Path | None = None, z_threshold: float = 3.0,
          module_file: Path | None = None):
    anndata_root = anndata_root if anndata_root is not None else output_root.parent / 'anndata'
    expression_path, attributes, patient, sample, inventory, legacy = acquire(cache)
    raw, features = read_rsem(expression_path)
    expression_all, raw_kept, feature_audit, stats = clean_expression(raw, features)
    audit = select_samples(expression_all, z_threshold)
    clinical, samples, audit, dup_stats = match_clinical(audit, patient, sample, inventory, legacy)
    order = samples.sample_id.tolist()
    expression = expression_all[order]
    rsem = raw_kept[order]
    genes = feature_audit.loc[feature_audit.selected, ['feature_id', 'gene_id', 'gene_symbol', 'source_row']].reset_index(drop=True)
    genes['mean_log2'] = expression.mean(axis=1).round(DECIMALS).to_numpy()
    genes['pct_samples_expressed'] = ((expression > DETECTED_LOG2).mean(axis=1) * 100).round(1).to_numpy()
    coverage = module_coverage(module_file, genes, feature_audit) if module_file is not None else None
    checks = validate(expression, clinical, samples, genes, feature_audit, audit)
    checks.update(dup_stats)
    checks['negative_rsem_values_clipped'] = stats['negative_rsem_values_clipped']
    metadata = build_metadata(cache, attributes, z_threshold, stats, module_file, clinical, samples, genes,
                              feature_audit, audit)

    output_root.mkdir(parents=True, exist_ok=True)
    anndata_root.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.tcga-', dir=output_root))
    anndata_stage = Path(tempfile.mkdtemp(prefix='.tcga-', dir=anndata_root))
    try:
        # Shortest round-trip float repr keeps the rounded log2 values compact and exact on reload.
        write_csv(expression, stage / 'expression.csv.gz', index=True, float_format=None)
        for frame, name in [(clinical, 'clinical.csv'), (samples, 'samples.csv'), (genes, 'genes.csv'),
                            (feature_audit, 'feature_audit.csv'), (audit, 'sample_audit.csv')]:
            write_csv(frame, stage / name, float_format=None)
        if coverage is not None:
            write_csv(coverage, stage / 'module_coverage.csv', float_format=None)
        reloaded = validate_export(stage, expression)
        if reloaded['counts'] != checks['counts']:
            raise ValueError('Export round-trip changed counts')
        checks['export_round_trip_valid'] = True
        write_anndata(expression, clinical, genes, metadata, anndata_stage / 'expression.h5ad', layers={'rsem': rsem})
        checks['anndata_round_trip_valid'] = True
        write_json(anndata_stage / 'metadata.json', metadata)
        (anndata_stage / 'checksums.sha256').write_text(''.join(
            f'{sha256(p)}  {p.name}\n' for p in sorted(anndata_stage.iterdir())))
        write_json(stage / 'metadata.json', metadata)
        write_json(stage / 'validation.json', checks)
        (stage / 'qc_summary.md').write_text(qc_summary(checks, stats, z_threshold))
        (stage / 'checksums.sha256').write_text(''.join(f'{sha256(p)}  {p.name}\n' for p in sorted(stage.iterdir())))
        for staged, root in [(stage, output_root), (anndata_stage, anndata_root)]:
            check_destination(staged, root / cache.root.name)
        publish(stage, output_root / cache.root.name)
        publish(anndata_stage, anndata_root / cache.root.name)
    finally:
        shutil.rmtree(stage, ignore_errors=True)
        shutil.rmtree(anndata_stage, ignore_errors=True)
    return checks


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--snapshot', default='tcga-prad-v1', help='Immutable input snapshot name; use a new name to refresh sources')
    parser.add_argument('--offline', action='store_true', help='Forbid HTTP; require all inputs and matching checksums')
    parser.add_argument('--output-root', type=Path, help='Alternate processed root, useful for offline verification')
    parser.add_argument('--anndata-root', type=Path, help='Alternate AnnData output root (default: data/anndata)')
    parser.add_argument('--z-threshold', type=float, default=3.0, help='QC outlier cutoff in SDs (default: 3)')
    parser.add_argument('--module-file', type=Path, help='Optional CSV/TSV of module genes (gene + module columns)')
    args = parser.parse_args(argv)
    if not args.snapshot or args.snapshot in {'.', '..'} or Path(args.snapshot).name != args.snapshot:
        parser.error('snapshot must be a single directory name')
    if args.z_threshold <= 0:
        parser.error('--z-threshold must be positive')
    cache = SourceCache(args.data_dir / 'raw' / args.snapshot, offline=args.offline,
                        accessions=[STUDY, LEGACY_STUDY], user_agent='TCGA-PRAD/1.0')
    output_root = args.output_root or args.data_dir / 'processed'
    anndata_root = args.anndata_root or args.data_dir / 'anndata'
    paths = [cache.root.resolve(), output_root.resolve(), anndata_root.resolve()]
    if any(a.is_relative_to(b) or b.is_relative_to(a) for i, a in enumerate(paths) for b in paths[i + 1:]):
        parser.error('raw, output and AnnData directories must not overlap')
    report = build(cache, output_root, anndata_root, args.z_threshold, args.module_file)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
