"""CLI for the Belfast / FASTMAN GSE116918 cohort (metastasis and BCR endpoints).

Reuses the shared GEO cohort helpers (feature selection, validation, CSV round trips, AnnData export and
output protection) without changing them. The shared builder supports one BCR endpoint, so this module adds
the metastasis endpoint and records it as the primary endpoint for module-scores.
"""
from __future__ import annotations

import argparse
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from . import belfast
from .anndata_export import write_anndata
from .expression import write_csv
from .geo_cohort_pipeline import finite_features, reload_table, validate as validate_bcr
from .sources import SourceCache, sha256, write_json
from .taylor_pipeline import check_destination, publish

FIELDS = ('event', 'months', 'eligibility_reason', 'eligible', 'time_unit')


def endpoint_counts(clinical, name):
    eligible = clinical[f'{name}_eligible']
    return {'eligible': int(eligible.sum()), 'events': int(clinical.loc[eligible, f'{name}_event'].eq(1).sum()),
            'censored': int(clinical.loc[eligible, f'{name}_event'].eq(0).sum()),
            'source_coded_events': int(clinical[f'{name}_event'].eq(1).sum()),
            'eligibility_reasons': clinical[f'{name}_eligibility_reason'].value_counts().to_dict()}


def validate(expression, clinical, samples, genes, mapping, audit):
    """Shared structural and BCR checks, plus the same re-decoding check for every other endpoint."""
    report = validate_bcr(expression, clinical, samples, genes, mapping, audit, belfast.outcome)
    for name in belfast.ENDPOINTS:
        if name == 'bcr':
            continue
        source = clinical.drop(columns=[f'{name}_{f}' for f in FIELDS])
        decoded = belfast.add_endpoint(source, name)
        for field in FIELDS:
            pd.testing.assert_series_equal(clinical[f'{name}_{field}'], decoded[f'{name}_{field}'], check_dtype=False)
    report['primary_endpoint'] = belfast.PRIMARY_ENDPOINT
    report['endpoints'] = {name: endpoint_counts(clinical, name) for name in belfast.ENDPOINTS}
    return report


def metadata_for(cache, clinical, samples, genes, processing, details):
    return {'schema_version': 1, 'study': belfast.STUDY, 'geo_series': belfast.SERIES, 'platform': belfast.PLATFORM,
        'snapshot': cache.root.name, 'source_manifest_sha256': sha256(cache.manifest_path),
        'primary_endpoint': belfast.PRIMARY_ENDPOINT,
        'expression': {'orientation': 'features_by_samples', 'scale': belfast.SCALE, 'source_processing': processing,
            'aggregation': 'none', 'additional_transformations': 'none',
            'missing_policy': 'Exclude whole feature when missing in any selected primary; never impute'},
        'endpoints': {name: {'definition': definition, 'unit': 'months', 'event_field': event, 'time_field': time,
                             'missing_policy': 'Retain nullable outcomes with eligibility reasons',
                             'zero_time_policy': 'Retain but exclude from survival eligibility; no pseudotime'}
                      for name, (event, time, definition) in belfast.ENDPOINTS.items()},
        'policies': details,
        'data_dictionary': {'expression.csv.gz': 'feature_id then sample IDs; deposited RMA log2 values',
            'genes.csv': 'One row per retained Almac probe: Entrez ID, symbol, Ensembl gene, probe orientation',
            'samples.csv': 'One selected biopsy per patient, in expression column order',
            'clinical.csv': 'Original GEO characteristics plus nullable metastasis_* and bcr_* outcomes',
            'feature_mapping.csv.gz': 'Every deposited probe, its annotation status and selection decision',
            'sample_audit.csv': 'Every GEO profile and its selection decision',
            'expression.h5ad': 'samples x features; clinical in obs, genes in var, JSON provenance in uns'},
        'output_schemas': {name: {c: str(f[c].dtype) for c in f}
                           for name, f in [('clinical.csv', clinical), ('samples.csv', samples), ('genes.csv', genes)]},
        'code_sha256': {p.name: sha256(p) for p in sorted(Path(__file__).parent.glob('*.py'))}}


def build(cache, output_root, interim_root, anndata_root):
    cache.verify_all()
    clinical, samples, audit, source, annotations, _, processing, details = belfast.prepare(cache)
    expression, genes, mapping = finite_features(source, annotations)
    report = validate(expression, clinical, samples, genes, mapping, audit)
    report['feature_dispositions'] = mapping.mapping_status.value_counts().to_dict()
    metadata = metadata_for(cache, clinical, samples, genes, processing, details)
    stages = []
    try:
        for root in (output_root, interim_root, anndata_root):
            root.mkdir(parents=True, exist_ok=True)
            stages.append(Path(tempfile.mkdtemp(prefix='.' + belfast.SNAPSHOT + '-', dir=root)))
        stage, interim_stage, ad_stage = stages
        write_csv(expression, stage / 'expression.csv.gz', index=True)
        tables = {'clinical.csv': clinical, 'samples.csv': samples, 'genes.csv': genes,
                  'feature_mapping.csv.gz': mapping, 'sample_audit.csv': audit}
        loaded = {}
        for name, frame in tables.items():
            write_csv(frame, stage / name)
            loaded[name] = reload_table(stage / name, frame)
        write_csv(source, interim_stage / 'feature_expression.csv.gz', index=True)
        saved = pd.read_csv(stage / 'expression.csv.gz', index_col=0, dtype={'feature_id': str}, float_precision='round_trip')
        if saved.index.tolist() != expression.index.tolist() or saved.columns.tolist() != expression.columns.tolist() \
                or not np.array_equal(saved.to_numpy(), expression.to_numpy()):
            raise ValueError('CSV round-trip changed expression values or labels')
        check = validate(saved, loaded['clinical.csv'], loaded['samples.csv'], loaded['genes.csv'],
                         loaded['feature_mapping.csv.gz'], loaded['sample_audit.csv'])
        if check['counts'] != report['counts'] or check['endpoints'] != report['endpoints']:
            raise ValueError('CSV round-trip changed counts')
        write_anndata(expression, clinical, genes, metadata, ad_stage / 'expression.h5ad')
        report.update(export_round_trip_valid=True, anndata_round_trip_valid=True)
        write_json(stage / 'validation.json', report)
        for destination in (stage, ad_stage):
            write_json(destination / 'metadata.json', metadata)
            (destination / 'checksums.sha256').write_text(
                ''.join(f'{sha256(p)}  {p.name}\n' for p in sorted(destination.iterdir())))
        for staged, root in zip(stages, (output_root, interim_root, anndata_root)):
            check_destination(staged, root / cache.root.name)
        for staged, root in zip(stages, (output_root, interim_root, anndata_root)):
            publish(staged, root / cache.root.name)
    finally:
        for stage in stages:
            shutil.rmtree(stage, ignore_errors=True)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=f'Acquire and export {belfast.STUDY} {belfast.SERIES}')
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--snapshot', default=belfast.SNAPSHOT)
    parser.add_argument('--offline', action='store_true')
    for name in ('output-root', 'interim-root', 'anndata-root'):
        parser.add_argument('--' + name, type=Path)
    args = parser.parse_args(argv)
    if not args.snapshot or args.snapshot in {'.', '..'} or Path(args.snapshot).name != args.snapshot:
        parser.error('snapshot must be a single directory name')
    roots = [args.output_root or args.data_dir / 'processed', args.interim_root or args.data_dir / 'interim',
             args.anndata_root or args.data_dir / 'anndata']
    paths = [p.resolve() for p in [args.data_dir / 'raw', *roots]]
    if any(a.is_relative_to(b) or b.is_relative_to(a) for i, a in enumerate(paths) for b in paths[i + 1:]):
        parser.error('raw and output directories must not overlap')
    cache = SourceCache(args.data_dir / 'raw' / args.snapshot, offline=args.offline,
                        accessions=[belfast.SERIES, belfast.PLATFORM], user_agent=belfast.SERIES + '/1.0')
    print(json.dumps(build(cache, *roots), indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
