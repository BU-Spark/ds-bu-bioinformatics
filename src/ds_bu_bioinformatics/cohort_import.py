"""Import explicitly curated research files; never claim to recover unavailable clinical data."""
from __future__ import annotations

import argparse
import csv
from functools import partial
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from .geo_cohort_pipeline import add_outcomes, build, decode_bcr, number
from .expression import open_text
from .sources import SourceCache, sha256

# Cohort identity is fixed here to prevent substituting a convenient, different study.
COHORTS = {
    'nakagawa': ('GSE10645', 'Nakagawa / Mayo'),
    'cpc-gene': ('GSE107299', 'CPC-GENE'),
    'cancermap': ('GSE94767', 'CancerMap'),
    'erho': ('GSE46691', 'Erho / Mayo'),
    'karnes': ('GSE62116', 'Karnes / Mayo'),
    'den': ('GSE72291', 'Den / Thomas Jefferson'),
    'klein': ('GSE62667', 'Klein / Cleveland'),
    'moscow': ('GSE229904', 'Moscow / EMB-RAS'),
    'dkfz-eopc': ('EGAD00001004791', 'DKFZ / EOPC'),
    'cpgea': ('PRJCA001124', 'CPGEA'),
}
FILES = ('expression', 'genes', 'samples', 'clinical', 'sample_audit')


def decode(record, prefix):
    unit = record.get(f'source__{prefix}_unit')
    if unit not in {'months', 'years'}:
        raise ValueError(f'{prefix}: every imported duration needs an explicit months/years unit')
    time = number(record.get(f'source__{prefix}_time'))
    if time is not None:
        time *= 12 if unit == 'years' else 1
        if not np.isfinite(time):
            raise ValueError('Imported duration overflows months')
    return decode_bcr(record.get(f'source__{prefix}_event'), time, {'0': 0, '1': 1})


def read_csv(path):
    # Keep identifiers, original missing tokens and leading zeroes intact.
    with open_text(path) as handle:
        header = next(csv.reader(handle), [])
    if not header or any(not value for value in header) or len(set(header)) != len(header):
        raise ValueError('CSV columns must be nonempty and unique; duplicate identifiers cannot be auto-renamed')
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def boolean(values):
    if not values.isin(['True', 'False', 'true', 'false', '1', '0']).all():
        raise ValueError('Audit selected must contain explicit true/false or 1/0')
    return values.isin(['True', 'true', '1'])


class BundleCohort:
    def __init__(self, key, manifest_path=None):
        self.key = key
        self.SERIES, self.STUDY = COHORTS[key]
        self.PLATFORM = 'curated-local-bundle'
        self.ACCESSIONS = [self.SERIES, 'curated-local-bundle-v1']
        self.SNAPSHOT = key + '-import-v1'
        self.manifest_path = manifest_path
        self.ENDPOINTS = {}

    def outcome(self, record):
        return next(iter(self.ENDPOINTS.values()))(record)

    def prepare(self, cache):
        path = cache.ingest('bundle/manifest.json', self.manifest_path, source='explicit curated cohort bundle')
        spec = json.loads(path.read_text())
        if spec.get('schema_version') != 1 or spec.get('study_accession') != self.SERIES:
            raise ValueError('Bundle schema/accession mismatch; cannot substitute another cohort')
        if spec.get('organism') != 'Homo sapiens':
            raise ValueError('Only explicitly annotated human primary cohorts are supported')
        for field in ('expression_scale', 'selection_policy', 'annotation_source', 'source_description'):
            if not isinstance(spec.get(field), str) or not spec[field].strip():
                raise ValueError(f'Bundle must document {field}')
        kind = spec.get('expression_kind')
        if kind not in {'raw_counts', 'tpm', 'fpkm', 'log_microarray', 'linear_microarray'}:
            raise ValueError('Specify expression_kind: raw_counts, tpm, fpkm, log_microarray or linear_microarray')
        platforms = spec.get('platforms')
        if not isinstance(platforms, list) or not platforms or any(not isinstance(p, str) or not p for p in platforms):
            raise ValueError('Bundle needs explicit platform names')
        definitions = spec.get('endpoints', {})
        primary = spec.get('primary_endpoint')
        if not isinstance(definitions, dict) or not definitions or primary not in definitions:
            raise ValueError('Bundle needs declared endpoints and a primary_endpoint')
        for prefix, definition in definitions.items():
            if not re.fullmatch(r'[a-z][a-z0-9]*', prefix):
                raise ValueError('Endpoint prefixes must be simple lower-case identifiers')
            if not isinstance(definition, dict) or not all(isinstance(definition.get(k), str) and definition[k].strip()
                    for k in ('definition', 'time_origin', 'censoring', 'sample_linkage')):
                raise ValueError(f'{prefix}: missing endpoint definition/origin/censoring/sample linkage')
        self.ENDPOINTS = {p: partial(decode, prefix=p) for p in [primary, *sorted(set(definitions) - {primary})]}
        self.PLATFORM = ';'.join(platforms)
        self.SCALE = spec['expression_scale']
        cached = {}
        for logical in FILES:
            entry = spec.get('files', {}).get(logical, {})
            relative = Path(entry.get('path', ''))
            digest = entry.get('sha256', '')
            if (not entry.get('path') or relative.is_absolute() or '..' in relative.parts
                    or not re.fullmatch(r'[0-9a-f]{64}', digest)
                    or not isinstance(entry.get('source'), str) or not entry['source'].strip()):
                raise ValueError(f'{logical}: specify a relative path, SHA-256 and original source citation')
            original = None
            if self.manifest_path is not None:
                original = self.manifest_path.parent / relative
                if not original.resolve().is_relative_to(self.manifest_path.parent.resolve()):
                    raise ValueError('Bundle input symlink escapes manifest directory')
            suffix = '.csv.gz' if logical == 'expression' else '.csv'
            cached[logical] = cache.ingest(f'bundle/{logical}{suffix}', original, source=entry['source'])
            if sha256(cached[logical]) != digest:
                raise ValueError(f'{logical}: supplied input checksum differs from bundle manifest')
        raw = read_csv(cached['expression'])
        if raw.empty or raw.columns[0] != 'feature_id' or raw.feature_id.eq('').any() or not raw.feature_id.is_unique:
            raise ValueError('Expression must have unique, nonempty original feature_id rows')
        source = raw.set_index('feature_id').apply(pd.to_numeric, errors='raise')
        # CSV blanks represent missing measurements; malformed nonempty tokens fail.
        genes = read_csv(cached['genes'])
        samples = read_csv(cached['samples'])
        clinical = read_csv(cached['clinical'])
        audit = read_csv(cached['sample_audit'])
        if not {'feature_id', 'gene_id', 'gene_symbol'}.issubset(genes) or not genes.feature_id.is_unique:
            raise ValueError('Gene table needs unique feature_id, gene_id and gene_symbol')
        if set(genes.feature_id) != set(source.index):
            raise ValueError('Gene annotations must cover exactly the source expression rows')
        genes = genes.set_index('feature_id').loc[source.index].reset_index()
        if not {'sample_id', 'patient_id', 'platform_id'}.issubset(samples) or samples.empty:
            raise ValueError('Samples need sample_id, patient_id and platform_id')
        if (not samples.sample_id.is_unique or not samples.patient_id.is_unique
                or samples.sample_id.eq('').any() or samples.patient_id.eq('').any()
                or not samples.platform_id.isin(platforms + [self.PLATFORM]).all()):
            raise ValueError('Missing/duplicate patient identities or undeclared platform')
        if not {'sample_id', 'selected', 'selection_reason'}.issubset(audit) or not audit.sample_id.is_unique:
            raise ValueError('Every source assay needs a unique sample-audit row')
        audit['selected'] = boolean(audit.selected)
        if (audit.selection_reason.eq('').any()
                or audit.loc[audit.selected, 'sample_id'].tolist() != samples.sample_id.tolist()
                or set(source.columns) != set(audit.sample_id)):
            raise ValueError('Source expression, complete sample audit and selected sample order disagree')
        if 'sample_id' not in clinical or not clinical.sample_id.is_unique or set(clinical.sample_id) != set(samples.sample_id):
            raise ValueError('Clinical and selected samples need exact one-to-one linkage')
        clinical = clinical.set_index('sample_id').loc[samples.sample_id].reset_index()
        for column in samples:
            if column in clinical and not clinical[column].equals(samples[column].reset_index(drop=True)):
                raise ValueError(f'Clinical/sample identity conflict: {column}')
            clinical[column] = samples[column].to_numpy()
        # Preserve imported original tokens before creating the shared endpoint columns.
        for prefix, decoder in self.ENDPOINTS.items():
            cols = [f'{prefix}_{suffix}' for suffix in ('event', 'time', 'unit')]
            if not set(cols).issubset(clinical):
                raise ValueError(f'{prefix}: provide event/time/unit columns, explicitly blank when unavailable')
            clinical = clinical.rename(columns={c: 'source__' + c for c in cols})
            clinical = add_outcomes(clinical, decoder, prefix)
            clinical[f'{prefix}_clinical_source'] = spec['source_description']
        # Never pass through unvalidated extra *_event columns to the shared checker.
        unexpected = [c for c in clinical if c.endswith('_event') and '__' not in c
                      and c[:-6] not in self.ENDPOINTS]
        if unexpected:
            raise ValueError(f'Undeclared endpoint columns: {unexpected}')
        source = source.loc[:, samples.sample_id].copy()
        library_sizes = None
        if kind != 'log_microarray' and source.lt(0).any().any():
            raise ValueError(f'{kind}: negative expression is invalid')
        if kind == 'raw_counts':
            values = source.to_numpy(dtype=float)
            if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
                raise ValueError('Raw counts must be complete finite nonnegative integers')
            totals = source.sum(axis=0)
            if not np.isfinite(totals).all() or not totals.gt(0).all():
                raise ValueError('Every raw-count sample needs a finite positive library total')
            # Freeze totals BEFORE annotation filtering so discarded features do
            # not silently change the CPM denominator used during module scoring.
            library_sizes = {str(k): float(v) for k, v in totals.items()}
        symbols = genes.gene_symbol
        valid = symbols.str.fullmatch(r'[^;/,\s]+') & ~symbols.str.upper().isin({'NA', 'N/A', '---', '-', 'NULL'})
        pairs = genes.loc[valid & genes.gene_id.ne(''), ['gene_id', 'gene_symbol']].drop_duplicates()
        bad_ids = set(pairs.loc[pairs.gene_id.duplicated(False), 'gene_id'])
        bad_symbols = set(pairs.loc[pairs.gene_symbol.duplicated(False), 'gene_symbol'])
        conflict = valid & (genes.gene_id.isin(bad_ids) | genes.gene_symbol.isin(bad_symbols))
        genes['mapping_status'] = np.select([conflict, valid],
            ['gene_symbol_identity_conflict', 'supplied_single_symbol'], default='invalid_or_ambiguous_annotation')
        genes['annotation_eligible'] = valid & ~conflict
        details = {k: spec[k] for k in ('selection_policy', 'annotation_source', 'source_description')}
        details.update(input_mode='curated_local_bundle', endpoints=definitions,
            expression_kind=kind, library_sizes=library_sizes,
            limitations=['Requires independently curated research inputs and correct source crosswalks; software checks cannot establish clinical validity',
                        'No direct parsing or independent verification of the original study files is claimed',
                        'No missing dates inferred; no raw sequencing/array preprocessing or probe recovery'],
            supplied_bundle_sha256=sha256(path))
        return clinical, samples, audit, source, genes, {}, [spec['source_description']], details

    def summarize(self, clinical, audit, mapping, extra):
        return {'input_mode': 'curated_local_bundle',
                'survival_eligibility': {p: int(clinical[f'{p}_eligible'].sum()) for p in self.ENDPOINTS}}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort', required=True, choices=sorted(COHORTS))
    parser.add_argument('--bundle', type=Path, help='First import: JSON manifest next to curated files; cached offline replay needs no originals')
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--snapshot')
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--clinical-file', type=Path)
    parser.add_argument('--clinical-provenance', type=Path)
    for name in ('output-root', 'interim-root', 'anndata-root'):
        parser.add_argument('--' + name, type=Path)
    args = parser.parse_args(argv)
    if (args.clinical_file is None) != (args.clinical_provenance is None):
        parser.error('supply both --clinical-file and --clinical-provenance')
    cohort = BundleCohort(args.cohort, args.bundle)
    snapshot = args.snapshot or cohort.SNAPSHOT
    if not snapshot or snapshot in {'.', '..'} or Path(snapshot).name != snapshot:
        parser.error('snapshot must be a single directory name')
    roots = [args.output_root or args.data_dir / 'processed', args.interim_root or args.data_dir / 'interim',
             args.anndata_root or args.data_dir / 'anndata']
    paths = [p.resolve() for p in [args.data_dir / 'raw', *roots]]
    if any(a.is_relative_to(b) or b.is_relative_to(a) for i, a in enumerate(paths) for b in paths[i + 1:]):
        parser.error('raw and output directories must not overlap')
    cache = SourceCache(args.data_dir / 'raw' / snapshot, offline=args.offline, accessions=cohort.ACCESSIONS)
    print(json.dumps(build(cache, *roots, cohort=cohort, clinical_file=args.clinical_file,
                           clinical_provenance=args.clinical_provenance), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
