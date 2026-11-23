"""Check every processed dataset against the shared output contract and summarize them in one table.

Each dataset has its own pipeline (taylor-gse21034, tcga-prad, cambridge-gse70768, ...). This command
does not rerun them. It reads data/processed/<snapshot>/ and data/anndata/<snapshot>/ and checks that
every dataset follows the same contract, so the outputs can be used together (for example for PCA
module scoring):

  * required files are present, and checksums.sha256 matches every file in each folder;
  * expression.csv.gz is feature_id x sample_id with unique labels and finite values;
  * samples.csv, clinical.csv and the selected rows of sample_audit.csv list the expression samples
    in the same order, with one sample per patient;
  * genes.csv lists the expression features in the same order;
  * every endpoint <e>_event column has <e>_months, <e>_eligible, <e>_eligibility_reason and
    <e>_time_unit; events are 0/1/missing; eligible rows have an event and a positive time;
  * validation.json reports a valid structure and both round trips, and metadata.json names the snapshot;
  * expression.h5ad has the same samples, features and values, and the same metadata in uns.

Counts in the summary are computed from the files, because validation.json count names differ
between datasets. Warnings flag differences that do not break the contract but matter for combined
analysis, such as expression that is not on a log scale.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .sources import sha256

PROCESSED_FILES = ('expression.csv.gz', 'samples.csv', 'clinical.csv', 'genes.csv', 'sample_audit.csv',
                   'metadata.json', 'validation.json', 'checksums.sha256')
FEATURE_AUDIT_FILES = ('feature_mapping.csv.gz', 'feature_audit.csv')  # one of these, depending on the dataset
ANNDATA_FILES = ('expression.h5ad', 'metadata.json', 'checksums.sha256')
REQUIRED_COLUMNS = {'samples.csv': ('sample_id', 'patient_id', 'platform_id'),
                    'clinical.csv': ('sample_id',),
                    'genes.csv': ('feature_id', 'gene_id', 'gene_symbol'),
                    'sample_audit.csv': ('sample_id', 'selected', 'selection_reason')}
ENDPOINT_SUFFIXES = ('months', 'eligible', 'eligibility_reason', 'time_unit')
TRUE_VALUES = {True, 'True', 'true', 1}
LOG_SCALE_MAX = 50.0  # log2 expression rarely exceeds ~25; larger values suggest an unlogged scale


class Check:
    """Collects failures for one dataset instead of stopping at the first one."""

    def __init__(self):
        self.errors: list[str] = []

    def require(self, condition, message: str) -> bool:
        if not condition:
            self.errors.append(message)
        return bool(condition)


def verify_checksums(folder: Path, check: Check, label: str) -> None:
    manifest = folder / 'checksums.sha256'
    listed = {}
    for line in manifest.read_text().splitlines():
        digest, _, name = line.partition('  ')
        listed[name] = digest
    present = {p.name for p in folder.iterdir() if p.is_file() and p.name != 'checksums.sha256'}
    check.require(set(listed) == present,
                  f'{label}: checksums.sha256 does not list exactly the files present '
                  f'(unlisted: {sorted(present - set(listed))}, absent: {sorted(set(listed) - present)})')
    bad = [name for name, digest in listed.items() if name in present and sha256(folder / name) != digest]
    check.require(not bad, f'{label}: checksum mismatch for {bad}')


def as_bool(series: pd.Series) -> pd.Series:
    return series.map(lambda value: value in TRUE_VALUES)


def check_endpoints(clinical: pd.DataFrame, check: Check) -> dict:
    prefixes = sorted(c[:-len('_event')] for c in clinical.columns if c.endswith('_event') and '__' not in c)
    check.require(prefixes, 'clinical.csv: no <endpoint>_event columns')
    summary = {}
    for prefix in prefixes:
        missing = [f'{prefix}_{s}' for s in ENDPOINT_SUFFIXES if f'{prefix}_{s}' not in clinical]
        if not check.require(not missing, f'clinical.csv: endpoint {prefix} lacks {missing}'):
            continue
        event = pd.to_numeric(clinical[f'{prefix}_event'], errors='coerce')
        months = pd.to_numeric(clinical[f'{prefix}_months'], errors='coerce')
        eligible = as_bool(clinical[f'{prefix}_eligible'])
        reason = clinical[f'{prefix}_eligibility_reason']
        check.require(clinical[f'{prefix}_event'].isna().eq(event.isna()).all() and event.dropna().isin([0, 1]).all(),
                      f'clinical.csv: {prefix}_event has values other than 0, 1 or missing')
        check.require(clinical[f'{prefix}_months'].isna().eq(months.isna()).all(),
                      f'clinical.csv: {prefix}_months has non-numeric values')
        check.require(np.isfinite(months.dropna()).all(),
                      f'clinical.csv: {prefix}_months has non-finite values')
        check.require(clinical[f'{prefix}_time_unit'].eq('months').fillna(False).all(),
                      f'clinical.csv: {prefix}_time_unit must be months for every row')
        check.require(eligible.eq(reason.eq('eligible')).all(),
                      f'clinical.csv: {prefix}_eligible disagrees with {prefix}_eligibility_reason')
        check.require(event[eligible].notna().all() and months[eligible].gt(0).all(),
                      f'clinical.csv: eligible {prefix} rows need an event and a positive time')
        check.require(reason.notna().all(), f'clinical.csv: {prefix}_eligibility_reason has missing values')
        summary[prefix] = {'eligible': int(eligible.sum()), 'events': int(event[eligible].eq(1).sum()),
                           'time_unit': ';'.join(sorted(clinical[f'{prefix}_time_unit'].dropna().astype(str).unique()))}
    return summary


def check_dataset(processed: Path, anndata_folder: Path, load_anndata: bool = True) -> dict:
    """Check one snapshot. Returns a summary dict with 'errors' (empty when the dataset passes)."""
    check = Check()
    result = {'dataset': processed.name, 'errors': check.errors, 'warnings': []}
    for name in PROCESSED_FILES:
        check.require((processed / name).is_file(), f'processed: missing {name}')
    check.require(any((processed / n).is_file() for n in FEATURE_AUDIT_FILES),
                  f'processed: missing a feature audit ({" or ".join(FEATURE_AUDIT_FILES)})')
    for name in ANNDATA_FILES:
        check.require((anndata_folder / name).is_file(), f'anndata: missing {name}')
    if check.errors:
        return result
    verify_checksums(processed, check, 'processed')
    verify_checksums(anndata_folder, check, 'anndata')

    expression = pd.read_csv(processed / 'expression.csv.gz', index_col=0, dtype={'feature_id': str}, float_precision='round_trip')
    tables = {name: pd.read_csv(processed / name, dtype={c: str for c in ('sample_id', 'patient_id', 'feature_id', 'gene_id')},
                                low_memory=False)
              for name in REQUIRED_COLUMNS}
    for name, columns in REQUIRED_COLUMNS.items():
        missing = [c for c in columns if c not in tables[name]]
        check.require(not missing, f'{name}: missing columns {missing}')
    if check.errors:
        return result
    samples, clinical, genes, audit = (tables[n] for n in ('samples.csv', 'clinical.csv', 'genes.csv', 'sample_audit.csv'))
    columns = expression.columns.astype(str).tolist()
    features = expression.index.astype(str).tolist()

    check.require(expression.index.name == 'feature_id', 'expression.csv.gz: index is not named feature_id')
    check.require(expression.index.is_unique and expression.columns.is_unique, 'expression.csv.gz: duplicate labels')
    values = expression.to_numpy()
    check.require(np.issubdtype(values.dtype, np.number) and np.isfinite(values).all(),
                  'expression.csv.gz: non-numeric or non-finite values')
    check.require(samples.sample_id.tolist() == columns, 'samples.csv: sample order differs from expression columns')
    check.require(clinical.sample_id.tolist() == columns, 'clinical.csv: sample order differs from expression columns')
    check.require(genes.feature_id.tolist() == features, 'genes.csv: feature order differs from expression rows')
    check.require(samples.patient_id.notna().all() and samples.patient_id.is_unique,
                  'samples.csv: missing or repeated patient_id (expected one sample per patient)')
    selected = audit.loc[as_bool(audit.selected), 'sample_id'].tolist()
    check.require(audit.sample_id.is_unique and sorted(selected) == sorted(columns),
                  'sample_audit.csv: selected samples differ from expression columns')
    endpoints = check_endpoints(clinical, check)

    validation = json.loads((processed / 'validation.json').read_text())
    for key in ('structurally_valid', 'export_round_trip_valid', 'anndata_round_trip_valid'):
        check.require(validation.get(key) is True, f'validation.json: {key} is not true')
    metadata = json.loads((processed / 'metadata.json').read_text())
    check.require(metadata.get('snapshot') == processed.name, 'metadata.json: snapshot does not match the folder name')
    check.require(json.loads((anndata_folder / 'metadata.json').read_text()) == metadata,
                  'anndata/metadata.json differs from processed/metadata.json')

    if load_anndata:
        import anndata as ad
        adata = ad.read_h5ad(anndata_folder / 'expression.h5ad')
        check.require(adata.obs_names.tolist() == columns, 'expression.h5ad: obs_names differ from expression columns')
        check.require(adata.var_names.tolist() == features, 'expression.h5ad: var_names differ from expression rows')
        x = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
        check.require(x.shape == values.T.shape and np.array_equal(x, values.T),
                      'expression.h5ad: X differs from expression.csv.gz')
        stored = adata.uns.get('metadata_json')
        check.require(stored is not None and json.loads(stored) == metadata,
                      "expression.h5ad: uns['metadata_json'] differs from metadata.json")
        for axis, filename, key in [('obs', 'clinical.csv', 'sample_id'), ('var', 'genes.csv', 'feature_id')]:
            original = pd.read_csv(processed / filename, dtype=str, keep_default_na=False).set_index(key)
            saved = getattr(adata, axis)
            check.require(original.columns.tolist() == saved.columns.tolist(),
                          f'expression.h5ad: {axis} annotation columns differ from {filename}')
            for column in original.columns.intersection(saved.columns):
                a, b = original[column].reset_index(drop=True), saved[column].reset_index(drop=True)
                if pd.api.types.is_numeric_dtype(b.dtype) and not pd.api.types.is_bool_dtype(b.dtype):
                    try:
                        a = pd.to_numeric(a.replace('', np.nan), errors='raise')
                        same = np.allclose(a.to_numpy(dtype=float), b.to_numpy(dtype=float, na_value=np.nan),
                                           rtol=1e-12, atol=0, equal_nan=True)
                    except (ValueError, TypeError):
                        same = False
                else:
                    same = a.astype('string').fillna('').equals(b.astype('string').fillna(''))
                check.require(same, f'expression.h5ad: {axis}[{column}] differs from {filename}')

    result.update({
        'samples': len(columns), 'patients': int(samples.patient_id.nunique()), 'feature_rows': len(features),
        'unique_genes': int(genes.gene_id.dropna().nunique()),
        'rows_without_gene_id': int(genes.gene_id.isna().sum()),
        'platform': ';'.join(sorted(samples.platform_id.dropna().astype(str).unique())),
        'expression_range': [float(values.min()), float(np.median(values)), float(values.max())],
        'endpoints': endpoints})
    if values.size and float(values.max()) > LOG_SCALE_MAX:
        result['warnings'].append(f'expression max is {float(values.max()):,.0f}, so values do not look log-transformed; '
                                  'other datasets are log2, so transform before PCA or combined analysis')
    return result


def discover(data_dir: Path, snapshots=None) -> list[str]:
    root = data_dir / 'processed'
    names = snapshots or sorted(p.name for p in root.iterdir() if p.is_dir() and not p.name.startswith('.'))
    return list(names)


def markdown(results: list[dict]) -> str:
    lines = ['| Dataset | Status | Samples | Feature rows | Unique genes | Platform | Expression min / median / max | Endpoints (eligible, events) |',
             '|---|---|---:|---:|---:|---|---|---|']
    for r in results:
        status = 'PASS' if not r['errors'] else f"FAIL ({len(r['errors'])})"
        if 'samples' not in r:
            lines.append(f"| {r['dataset']} | {status} | | | | | | |")
            continue
        low, mid, high = r['expression_range']
        endpoints = '; '.join(f"{p.upper()} {e['eligible']}, {e['events']}" for p, e in r['endpoints'].items())
        lines.append(f"| {r['dataset']} | {status} | {r['samples']} | {r['feature_rows']:,} | {r['unique_genes']:,} | "
                     f"{r['platform']} | {low:.2f} / {mid:.2f} / {high:.2f} | {endpoints} |")
    failures = [(r['dataset'], e) for r in results for e in r['errors']]
    if failures:
        lines += ['', 'Failures:'] + [f'- {dataset}: {error}' for dataset, error in failures]
    warnings = [(r['dataset'], w) for r in results for w in r.get('warnings', [])]
    if warnings:
        lines += ['', 'Warnings:'] + [f'- {dataset}: {warning}' for dataset, warning in warnings]
    return '\n'.join(lines) + '\n'


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--snapshot', action='append', help='Check only this snapshot (repeatable); default: all')
    parser.add_argument('--skip-anndata', action='store_true', help='Do not load expression.h5ad files')
    parser.add_argument('--json', type=Path, help='Also write the full results as JSON')
    args = parser.parse_args(argv)
    if not (args.data_dir / 'processed').is_dir():
        parser.error(f'{args.data_dir / "processed"} does not exist')
    results = []
    for name in discover(args.data_dir, args.snapshot):
        try:
            results.append(check_dataset(args.data_dir / 'processed' / name, args.data_dir / 'anndata' / name,
                                         load_anndata=not args.skip_anndata))
        except Exception as error:  # report unreadable files as a failure for that dataset, then continue
            results.append({'dataset': name, 'errors': [f'could not be read: {type(error).__name__}: {error}'],
                            'warnings': []})
    print(markdown(results), end='')
    if args.json:
        args.json.write_text(json.dumps(results, indent=2, sort_keys=True) + '\n')
    return 1 if any(r['errors'] for r in results) else 0


if __name__ == '__main__':
    raise SystemExit(main())
