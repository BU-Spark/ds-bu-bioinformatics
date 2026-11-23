"""Shared finite-feature CSV/H5AD export with explicit cohort endpoint adapters."""
from __future__ import annotations
from .provenance import code_hashes

import argparse
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .anndata_export import write_anndata
from .expression import write_csv
from .sources import SourceCache, sha256, write_json
from .taylor_pipeline import check_destination, publish

MISSING = {'', 'NA', 'N/A', 'UNKNOWN', '[NOT AVAILABLE]'}


def number(value):
    if pd.isna(value) or str(value).strip().upper() in MISSING:
        return None
    result = float(value)
    if not np.isfinite(result) or result < 0:
        raise ValueError(f'Invalid nonnegative clinical value: {value!r}')
    return result


def decode_bcr(status, duration, codes):
    months = number(duration)
    if pd.isna(status) or str(status).strip().upper() in MISSING:
        return None, months, 'missing_status'
    if status not in codes:
        raise ValueError(f'Unknown BCR status: {status!r}')
    event = codes[status]
    reason = 'missing_time' if months is None else 'zero_time_requires_review' if months == 0 else 'eligible'
    return event, months, reason


def add_outcomes(clinical, outcome, prefix='bcr'):
    clinical = clinical.copy()
    rows = [outcome(r) for r in clinical.to_dict('records')]
    clinical[f'{prefix}_event'] = pd.array([r[0] for r in rows], dtype='Int8')
    clinical[f'{prefix}_months'] = pd.array([r[1] for r in rows], dtype='Float64')
    clinical[f'{prefix}_eligibility_reason'] = [r[2] for r in rows]
    clinical[f'{prefix}_eligible'] = clinical[f'{prefix}_eligibility_reason'].eq('eligible')
    clinical[f'{prefix}_time_unit'] = 'months'
    return clinical


def finite_features(source, annotations):
    if source.index.tolist() != annotations.feature_id.tolist() or not source.index.is_unique:
        raise ValueError('Feature annotations are not aligned or unique')
    mapping = annotations.copy()
    mapping['missing_selected_values'] = source.isna().sum(axis=1).to_numpy()
    mapping['expression_status'] = np.where(mapping.missing_selected_values.gt(0), 'missing_selected_expression', 'complete')
    if np.isinf(source.to_numpy()).any():
        raise ValueError('Infinite source expression')
    mapping['selected'] = mapping.annotation_eligible & mapping.expression_status.eq('complete')
    genes = mapping.loc[mapping.selected, annotations.columns].drop(columns=['annotation_eligible', 'mapping_status']).reset_index(drop=True)
    expression = source.loc[genes.feature_id].copy()
    if expression.empty:
        raise ValueError('No finite eligible expression features')
    return expression, genes, mapping


def validate(expression, clinical, samples, genes, mapping, audit, outcome, *, endpoints=None):
    if expression.empty or not np.isfinite(expression.to_numpy()).all():
        raise ValueError('Empty or nonfinite final expression')
    if not expression.index.is_unique or not expression.columns.is_unique:
        raise ValueError('Duplicate expression labels')
    if expression.columns.tolist() != samples.sample_id.tolist() or expression.columns.tolist() != clinical.sample_id.tolist():
        raise ValueError('Expression/sample/clinical order differs')
    if not samples.equals(clinical[samples.columns]) or samples.patient_id.isna().any() or not samples.patient_id.is_unique:
        raise ValueError('Missing/duplicate patients or clinical/sample identity mismatch')
    if expression.index.tolist() != genes.feature_id.tolist() or not genes.feature_id.is_unique:
        raise ValueError('Expression/feature annotation order differs')
    if not mapping.feature_id.is_unique or mapping.loc[mapping.selected, 'feature_id'].tolist() != expression.index.tolist():
        raise ValueError('Feature audit does not match selection')
    expected_genes = mapping.loc[mapping.selected, genes.columns].reset_index(drop=True)
    pd.testing.assert_frame_equal(genes, expected_genes, check_dtype=False)
    if not audit.sample_id.is_unique or audit.loc[audit.selected, 'sample_id'].tolist() != samples.sample_id.tolist():
        raise ValueError('Sample audit does not match selection')
    endpoints = endpoints or {'bcr': outcome}
    counts, reasons = {}, {}
    for prefix, decode in endpoints.items():
        fields = [f'{prefix}_{suffix}' for suffix in ['event', 'months', 'eligibility_reason', 'eligible', 'time_unit']]
        decoded = add_outcomes(clinical.drop(columns=fields), decode, prefix)
        for field in fields:
            pd.testing.assert_series_equal(clinical[field], decoded[field], check_dtype=False)
        eligible = clinical[f'{prefix}_eligible']
        event = clinical[f'{prefix}_event']
        counts.update({f'{prefix}_events': int(event[eligible].eq(1).sum()),
            f'{prefix}_censored': int(event[eligible].eq(0).sum()),
            f'eligible_{prefix}_profiles': int(eligible.sum()),
            f'source_coded_{prefix}_events': int(event.eq(1).sum()),
            f'source_coded_{prefix}_non_events': int(event.eq(0).sum())})
        reasons[f'{prefix}_eligibility_reasons'] = clinical[f'{prefix}_eligibility_reason'].value_counts().to_dict()
    primary = next(iter(endpoints))
    return {'structurally_valid': True, 'counts': {
        'primary_profiles': len(samples), 'retained_feature_rows': len(genes),
        'unique_genes': int(genes.gene_id.nunique()),
        **counts, 'censored': counts[f'{primary}_censored'],
        'source_coded_non_events': counts[f'source_coded_{primary}_non_events']},
        'geo_profiles': len(audit), 'source_features': len(mapping),
        'sample_dispositions': audit.selection_reason.value_counts().to_dict(),
        'feature_dispositions': mapping.mapping_status.value_counts().to_dict(),
        'features_with_missing_selected_expression': int(mapping.missing_selected_values.gt(0).sum()),
        'missing_selected_expression_cells': int(mapping.missing_selected_values.sum()),
        'otherwise_eligible_features_excluded_for_missing_expression': int((mapping.annotation_eligible & mapping.missing_selected_values.gt(0)).sum()),
        **reasons,
        'clinical_missingness': clinical.isna().sum().to_dict(),
        'expression_range': [float(expression.to_numpy().min()), float(expression.to_numpy().max())]}


def reload_table(path, expected):
    # Preserve identifiers and original NA/N/A text as strings. Only an empty
    # serialized field represents null; clinical numeric columns remain numeric.
    dtypes = {c: 'string' for c in expected if pd.api.types.is_string_dtype(expected[c]) or expected[c].dtype == object}
    for c in expected:
        if isinstance(expected[c].dtype, pd.Int8Dtype | pd.Int64Dtype):
            dtypes[c] = str(expected[c].dtype)
    saved = pd.read_csv(path, dtype=dtypes, keep_default_na=False, na_values=[''], float_precision='round_trip')
    for c in expected:
        if isinstance(expected[c].dtype, pd.Float64Dtype):
            saved[c] = saved[c].astype('Float64')
        elif pd.api.types.is_bool_dtype(expected[c]):
            if len(saved) == 0:
                saved[c] = saved[c].astype(bool)
            if saved[c].dtype != bool:
                raise ValueError('CSV changed boolean field ' + c)
        if c in dtypes and dtypes[c] == 'string':
            pd.testing.assert_series_equal(saved[c].fillna(''), expected[c].astype('string').fillna(''), check_dtype=False)
        else:
            pd.testing.assert_series_equal(saved[c], expected[c], check_dtype=False)
    return saved


def build(cache, output_root, interim_root, anndata_root, *, cohort, clinical_file=None, clinical_provenance=None):
    with cache.transaction():
        return _build(cache, output_root, interim_root, anndata_root, cohort=cohort,
                      clinical_file=clinical_file, clinical_provenance=clinical_provenance)


def _build(cache, output_root, interim_root, anndata_root, *, cohort, clinical_file=None, clinical_provenance=None):
    cache.verify_all()
    prepared = cohort.prepare(cache)
    clinical, samples, audit, source, annotations, extra, processing, details = prepared
    expression, genes, mapping = finite_features(source, annotations)
    endpoints = getattr(cohort, 'ENDPOINTS', {'bcr': cohort.outcome})
    linkage_metadata = None
    from . import clinical_linkage
    if (clinical_file is None) != (clinical_provenance is None):
        raise ValueError('Supply both clinical_file and clinical_provenance')
    if clinical_file is not None or clinical_linkage.RECORD_NAME in cache.manifest['sources']:
        record_path = cache.ingest(clinical_linkage.RECORD_NAME, clinical_file, source='explicit clinical supplement')
        provenance_path = cache.ingest(clinical_linkage.PROVENANCE_NAME, clinical_provenance, source='explicit clinical supplement provenance')
        clinical, endpoints, linked_records, linkage_metadata = clinical_linkage.apply(
            clinical, endpoints, record_path, provenance_path, accession=cohort.SERIES,
            source_definitions=details.get('endpoints', {}))
        extra = {**extra, 'clinical_linkage_audit.csv': linked_records}
    report = validate(expression, clinical, samples, genes, mapping, audit, cohort.outcome, endpoints=endpoints)
    report.update(cohort.summarize(clinical, audit, mapping, extra))
    from .analysis_policy import policy
    analysis_policy = policy(cohort.SERIES)
    report.pop('survival_ready', None)
    report['complete_clinical_records'] = {p: int(clinical[f'{p}_eligible'].sum()) for p in endpoints}
    report['analysis_policy'] = analysis_policy
    metadata = {'schema_version': 1, 'study': cohort.STUDY, 'geo_series': cohort.SERIES, 'platform': cohort.PLATFORM,
        'snapshot': cache.root.name, 'source_manifest_sha256': sha256(cache.manifest_path),
        'expression': {'orientation': 'features_by_samples', 'scale': cohort.SCALE,
            'source_processing': processing, 'aggregation': 'none',
            'additional_transformations': details.get('normalization', 'none'),
            'missing_policy': 'Exclude whole feature when missing in any selected primary; never impute'},
        'bcr': {'endpoint': 'source biochemical recurrence', 'unit': 'months',
            'source': details.get('bcr_source', 'GEO sample characteristics'), 'missing_policy': 'Retain nullable outcomes with eligibility reasons',
            'zero_time_policy': 'Retain but exclude from survival eligibility; no pseudotime',
            'event_counts_policy': 'bcr_events/censored count eligible observations; source_coded_* report all known statuses'},
        'policies': details, 'analysis_policy': analysis_policy,
        'data_dictionary': {'expression.csv.gz': 'feature_id then sample IDs; source numeric values',
            'genes.csv': 'One row per expression feature, with original source identifiers and annotations',
            'samples.csv': 'One selected assay per patient, in expression column order',
            'clinical.csv': 'Original GEO/source clinical fields plus nullable bcr_* outcomes and eligibility',
            'feature_mapping.csv.gz': 'Every source expression row, annotation/missingness and selection decision',
            'sample_audit.csv': 'Every GEO profile, source identity and sample selection decision',
            'expression.h5ad': 'samples x features; clinical in obs, genes in var, JSON provenance in uns'},
        'output_schemas': {name: {c: str(f[c].dtype) for c in f} for name, f in [('clinical.csv', clinical), ('samples.csv', samples), ('genes.csv', genes), *extra.items()]},
        'code_sha256': code_hashes(__name__, cohort)}
    if hasattr(cohort, 'ENDPOINTS'):
        metadata.pop('bcr')
        metadata['primary_endpoint'] = next(iter(endpoints))
        metadata['endpoints'] = details['endpoints']
        metadata['data_dictionary']['clinical.csv'] = 'Original source fields plus nullable endpoint event/months/eligibility fields; unavailable times are never inferred'
    if linkage_metadata is not None:
        metadata['clinical_linkage'] = linkage_metadata
        # Keep the deposited definition authoritative; row-specific supplemental
        # definitions/provenance stay separate and never relabel unmatched rows.
        metadata['clinical_linkage']['applies_to'] = 'Matched rows only; see linked__* fields and event_source/duration_source'
        metadata.setdefault('primary_endpoint', next(iter(endpoints)))
    if details.get('expression_kind'):
        metadata['expression']['kind'] = details['expression_kind']
        if details.get('library_sizes') is not None:
            metadata['expression']['library_sizes'] = details['library_sizes']
    metadata['expression']['downstream_preprocessing'] = (
        'Exports retain original source values. Consumers must normalize for their analysis using '
        'the declared expression kind and full-source library totals for raw counts. '
        'This exporter does not normalize for analysis or fit survival models.')
    stages = []
    try:
        for root in (output_root, interim_root, anndata_root):
            root.mkdir(parents=True, exist_ok=True)
            stages.append(Path(tempfile.mkdtemp(prefix='.' + cohort.SNAPSHOT + '-', dir=root)))
        stage, interim_stage, ad_stage = stages
        write_csv(expression, stage / 'expression.csv.gz', index=True)
        tables = {'clinical.csv': clinical, 'samples.csv': samples, 'genes.csv': genes,
                  'feature_mapping.csv.gz': mapping, 'sample_audit.csv': audit, **extra}
        loaded = {}
        for name, frame in tables.items():
            write_csv(frame, stage / name)
            loaded[name] = reload_table(stage / name, frame)
        write_csv(source, interim_stage / 'feature_expression.csv.gz', index=True)
        (interim_stage / 'checksums.sha256').write_text(
            f'{sha256(interim_stage / "feature_expression.csv.gz")}  feature_expression.csv.gz\n')
        saved_expression = pd.read_csv(stage / 'expression.csv.gz', index_col=0, dtype={'feature_id': str}, float_precision='round_trip')
        if saved_expression.index.tolist() != expression.index.tolist() or saved_expression.columns.tolist() != expression.columns.tolist() or not np.array_equal(saved_expression.to_numpy(), expression.to_numpy()):
            raise ValueError('CSV round-trip changed expression values or labels')
        check = validate(saved_expression, loaded['clinical.csv'], loaded['samples.csv'], loaded['genes.csv'], loaded['feature_mapping.csv.gz'], loaded['sample_audit.csv'], cohort.outcome, endpoints=endpoints)
        if check['counts'] != report['counts']:
            raise ValueError('CSV round-trip changed counts')
        write_anndata(expression, clinical, genes, metadata, ad_stage / 'expression.h5ad')
        report.update(export_round_trip_valid=True, anndata_round_trip_valid=True)
        write_json(stage / 'validation.json', report)
        for destination in (stage, ad_stage):
            write_json(destination / 'metadata.json', metadata)
            (destination / 'checksums.sha256').write_text(''.join(f'{sha256(p)}  {p.name}\n' for p in sorted(destination.iterdir())))
        for staged, root in zip(stages, (output_root, interim_root, anndata_root)):
            check_destination(staged, root / cache.root.name)
        for staged, root in zip(stages, (output_root, interim_root, anndata_root)):
            publish(staged, root / cache.root.name)
    finally:
        for stage in stages:
            shutil.rmtree(stage, ignore_errors=True)
    return report


def main(cohort, argv=None):
    parser = argparse.ArgumentParser(description=f'Acquire and export {cohort.STUDY} {cohort.SERIES}')
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--snapshot', default=cohort.SNAPSHOT)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--clinical-file', type=Path, help='Verified sample-linked endpoint CSV; must also supply --clinical-provenance')
    parser.add_argument('--clinical-provenance', type=Path, help='JSON source, endpoint definitions, time origins, censoring and sample linkage')
    for name in ('output-root', 'interim-root', 'anndata-root'):
        parser.add_argument('--' + name, type=Path)
    args = parser.parse_args(argv)
    if (args.clinical_file is None) != (args.clinical_provenance is None):
        parser.error('supply both --clinical-file and --clinical-provenance')
    if not args.snapshot or args.snapshot in {'.', '..'} or Path(args.snapshot).name != args.snapshot:
        parser.error('snapshot must be a single directory name')
    roots = [args.output_root or args.data_dir/'processed', args.interim_root or args.data_dir/'interim', args.anndata_root or args.data_dir/'anndata']
    paths = [p.resolve() for p in [args.data_dir/'raw', *roots]]
    if any(a.is_relative_to(b) or b.is_relative_to(a) for i, a in enumerate(paths) for b in paths[i+1:]):
        parser.error('raw and output directories must not overlap')
    cache = SourceCache(args.data_dir/'raw'/args.snapshot, offline=args.offline,
        accessions=getattr(cohort, 'ACCESSIONS', [cohort.SERIES, cohort.PLATFORM]), user_agent=cohort.SERIES+'/1.0')
    print(json.dumps(build(cache, *roots, cohort=cohort,
        clinical_file=args.clinical_file, clinical_provenance=args.clinical_provenance), indent=2, sort_keys=True))
