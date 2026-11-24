"""Study-bound DKFZ/CPGEA workflows for authorized processed RNA-seq source files.

These parse expression and original sample/clinical tables directly using an explicit
column map, perform tumor selection/annotation, and generate every shared export.
They do not download controlled FASTQs or require a prebuilt five-file cohort bundle.
"""
import argparse
from functools import partial
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .gene_identifiers import annotate
from .geo_cohort_pipeline import add_outcomes, build, decode_bcr, number
from .geo_io import align_columns, read_tsv
from .moscow import identifier_scope
from .sources import SourceCache

STUDIES = {'dkfz-eopc': ('EGAD00001004791', 'DKFZ / EOPC'), 'cpgea': ('PRJCA001124', 'CPGEA')}


def source_field(column):
    """Reversible source-header escaping for HDF5 keys; percent escapes first."""
    return 'source__' + column.replace('%', '%25').replace('/', '%2F')


def source_outcome(record, prefix, definition):
    status_column = definition.get('event_column')
    time_column = definition.get('time_column')
    status = record.get(source_field(status_column)) if status_column else None
    time = number(record.get(source_field(time_column))) if time_column else None
    unit = definition.get('time_unit')
    if unit not in {'months', 'years'}:
        raise ValueError('Authorized endpoint needs explicit months/years units')
    if time is not None:
        time *= 12 if unit == 'years' else 1
        if not np.isfinite(time):
            raise ValueError('Clinical duration conversion overflows months')
    return decode_bcr(status, time, definition.get('event_codes', {}))


class AuthorizedRNASeq:
    def __init__(self, key, *, expression_file=None, sample_file=None, source_config=None):
        self.key = key
        self.SERIES, self.STUDY = STUDIES[key]
        self.PLATFORM = 'Illumina RNA-seq'
        self.ACCESSIONS = [self.SERIES, 'authorized-processed-rnaseq-v1']
        self.SNAPSHOT = key + '-v1'
        self.expression_file, self.sample_file, self.source_config = expression_file, sample_file, source_config
        self.ENDPOINTS = {}

    def outcome(self, record):
        return next(iter(self.ENDPOINTS.values()))(record)

    def cached_table(self, cache, logical, original):
        existing = [name for name in (f'authorized/{logical}.tsv', f'authorized/{logical}.tsv.gz') if name in cache.manifest['sources']]
        if len(existing) > 1:
            raise ValueError('Multiple authorized source versions in snapshot')
        name = existing[0] if existing else f'authorized/{logical}.tsv' + ('.gz' if original and original.suffix == '.gz' else '')
        return cache.ingest(name, original, source=self.SERIES + ' authorized ' + logical + ' source')

    def prepare(self, cache):
        config_path = cache.ingest('authorized/source-config.json', self.source_config,
                                  source=self.SERIES + ' source schema/provenance')
        config = json.loads(config_path.read_text())
        if config.get('study_accession') != self.SERIES:
            raise ValueError('Wrong study accession; a different DKFZ/Asian cohort cannot be substituted')
        for key in ('source', 'expression_scale', 'sample_linkage', 'specimen_definition'):
            if not isinstance(config.get(key), str) or not config[key].strip():
                raise ValueError('Source config must document ' + key)
        kind = config.get('expression_kind')
        if kind not in {'raw_counts', 'tpm', 'fpkm'}:
            raise ValueError('Specify raw_counts, tpm or fpkm for authorized processed expression')
        definitions = config.get('endpoints', {})
        primary = config.get('primary_endpoint')
        if not definitions or primary not in definitions:
            raise ValueError('Declare source endpoints and primary_endpoint; unavailable outcomes stay explicitly missing')
        for prefix, definition in definitions.items():
            if not prefix.isalnum() or not prefix.islower():
                raise ValueError('Invalid endpoint prefix')
            if not all(isinstance(definition.get(k), str) and definition[k].strip() for k in ('definition', 'time_origin', 'censoring')):
                raise ValueError('Document endpoint meaning, time origin and censoring')
            codes = definition.get('event_codes', {})
            if not isinstance(codes, dict) or any(v not in (0, 1) for v in codes.values()):
                raise ValueError('Source event codes must map explicitly to 0/1')
            if definition.get('time_unit') not in {'months', 'years'}:
                raise ValueError('Endpoint time unit must be months/years')
        self.ENDPOINTS = {p: partial(source_outcome, prefix=p, definition=definitions[p])
                          for p in [primary, *sorted(set(definitions) - {primary})]}
        self.SCALE = config['expression_scale']
        expression = read_tsv(self.cached_table(cache, 'expression', self.expression_file))
        sample_path = self.cached_table(cache, 'samples', self.sample_file)
        # No numeric coercion of clinical IDs, leading zeros or original missing labels.
        source_meta = pd.read_csv(sample_path, sep='\t', dtype=str, keep_default_na=False)
        columns = config.get('sample_columns', {})
        required = {'sample_id', 'patient_id', 'expression_id', 'tissue'}
        if not required.issubset(columns) or not set(columns.values()).issubset(source_meta):
            raise ValueError('Map sample_id, patient_id, expression_id and tissue to actual source columns')
        if source_meta.empty or not source_meta[columns['sample_id']].is_unique or not source_meta[columns['expression_id']].is_unique:
            raise ValueError('Missing/duplicate source sample or expression identity')
        for name in ('sample_id', 'patient_id', 'expression_id'):
            if source_meta[columns[name]].eq('').any():
                raise ValueError('Empty authorized source identity: ' + name)
        for definition in definitions.values():
            for key in ('event_column', 'time_column'):
                if definition.get(key) is not None and definition[key] not in source_meta:
                    raise ValueError('Configured clinical source field absent: ' + definition[key])
        source_columns = {c: source_field(c) for c in source_meta}
        meta = source_meta.rename(columns=source_columns)
        for target, original in columns.items():
            meta[target] = source_meta[original]
        meta['geo_accession'] = self.SERIES + ':' + meta.sample_id
        meta['platform_id'] = self.PLATFORM
        meta['title'] = meta.sample_id
        expression, expression_linkage = align_columns(expression, meta, lambda r: [r['expression_id']])
        tumor_values, excluded_values = set(config.get('tumor_values', [])), set(config.get('excluded_tissue_values', []))
        if not tumor_values or tumor_values & excluded_values or not meta.tissue.isin(tumor_values | excluded_values).all():
            raise ValueError('Classify all actual source specimen values explicitly as tumor or excluded')
        primary_mask = meta.tissue.isin(tumor_values)
        present = meta.sample_id.isin(expression.columns)
        candidates = meta.loc[primary_mask & present]
        if candidates.empty:
            raise ValueError('No primary tumors link to authorized expression columns')
        for prefix, decoder in self.ENDPOINTS.items():
            for patient, group in candidates.groupby('patient_id', sort=False):
                decoded = [decoder(r) for r in group.to_dict('records')]
                events = {r[0] for r in decoded if r[0] is not None}
                times = {r[1] for r in decoded if r[1] is not None}
                if len(events) > 1 or len(times) > 1:
                    raise ValueError(f'Conflicting {prefix} clinical records across tumor areas for patient {patient}')
        policy = config.get('replicate_policy', 'error')
        if policy not in {'error', 'lowest_sample_id'}:
            raise ValueError('Unsupported replicate policy')
        if candidates.patient_id.duplicated().any() and policy == 'error':
            raise ValueError('Multiple tumor areas/assays per patient; specify an explicit independent selection policy')
        winners = set(candidates.sort_values('sample_id').drop_duplicates('patient_id').sample_id)
        meta['selected'] = meta.sample_id.isin(winners)
        meta['selection_reason'] = np.select([~primary_mask, primary_mask & ~present, meta.selected],
            ['excluded_source_specimen', 'tumor_absent_from_expression', 'selected_primary_tumor'], default='other_region_same_patient')
        clinical = meta.loc[meta.selected].reset_index(drop=True)
        for prefix, decoder in self.ENDPOINTS.items():
            clinical = add_outcomes(clinical, decoder, prefix)
            clinical[f'{prefix}_clinical_source'] = config['source']
        samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title', 'expression_id']].copy()
        source = expression.loc[:, samples.sample_id]
        if source.lt(0).any().any():
            raise ValueError('Negative RNA-seq expression')
        totals = None
        if kind == 'raw_counts':
            values = source.to_numpy(dtype=float)
            if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
                raise ValueError('Raw counts require complete finite integer measurements')
            sizes = source.sum(axis=0)
            if not np.isfinite(sizes).all() or not sizes.gt(0).all():
                raise ValueError('Every selected sample needs a positive finite library total')
            totals = {str(k): float(v) for k, v in sizes.items()}
        annotations, gene_candidates = annotate(cache, source.index.tolist(), scope=identifier_scope(source.index))
        details = {'sample_selection': config['specimen_definition'] + '; replicate policy: ' + policy,
            'patient_identity': config['sample_linkage'], 'source_provenance': config,
            'source_column_export_names': source_columns,
            'annotation': 'Source gene IDs automatically validated through cached exact human identifier mapping; original rows preserved',
            'expression_kind': kind, 'library_sizes': totals,
            'endpoints': definitions,
            'limitations': ['Requires authorized processed expression and matching original sample/clinical tables; no controlled FASTQ acquisition or read alignment',
                'Full real-source files have not been validated in this environment; synthetic source-format tests only',
                'No literature-size assumptions, missing date inference or automatic merging of tumor areas',
                'Source column maps document interpretation; software cannot independently adjudicate clinical semantics']}
        return clinical, samples, meta, source, annotations, {'expression_linkage.csv': expression_linkage,
            'gene_mapping_candidates.csv.gz': gene_candidates}, [config['source']], details

    def summarize(self, clinical, audit, mapping, extra):
        return {'survival_ready': any(clinical[f'{p}_eligible'].any() for p in self.ENDPOINTS),
                'other_regions_excluded': int(audit.selection_reason.eq('other_region_same_patient').sum())}


def main(key, argv=None):
    parser = argparse.ArgumentParser(description=f'Export authorized processed RNA-seq sources for {STUDIES[key][1]}')
    parser.add_argument('--expression-file', type=Path, help='Original processed gene-by-assay TSV(.gz), first column gene ID')
    parser.add_argument('--sample-file', type=Path, help='Source sample/clinical TSV(.gz) with actual donor and specimen fields')
    parser.add_argument('--source-config', type=Path, help='Study-specific JSON column map, endpoint interpretation and source provenance')
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--snapshot', default=key + '-v1')
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--clinical-file', type=Path)
    parser.add_argument('--clinical-provenance', type=Path)
    for name in ('output-root', 'interim-root', 'anndata-root'):
        parser.add_argument('--' + name, type=Path)
    args = parser.parse_args(argv)
    supplied = [args.expression_file, args.sample_file, args.source_config]
    if any(v is not None for v in supplied) and not all(v is not None for v in supplied):
        parser.error('first acquisition requires --expression-file, --sample-file and --source-config together')
    if (args.clinical_file is None) != (args.clinical_provenance is None):
        parser.error('supply both clinical supplement flags')
    if args.snapshot in {'', '.', '..'} or Path(args.snapshot).name != args.snapshot:
        parser.error('snapshot must be a single directory name')
    roots = [args.output_root or args.data_dir / 'processed', args.interim_root or args.data_dir / 'interim',
             args.anndata_root or args.data_dir / 'anndata']
    paths = [p.resolve() for p in [args.data_dir / 'raw', *roots]]
    if any(a.is_relative_to(b) or b.is_relative_to(a) for i, a in enumerate(paths) for b in paths[i + 1:]):
        parser.error('raw and output roots must not overlap')
    cohort = AuthorizedRNASeq(key, expression_file=args.expression_file, sample_file=args.sample_file, source_config=args.source_config)
    cache = SourceCache(args.data_dir / 'raw' / args.snapshot, offline=args.offline, accessions=cohort.ACCESSIONS)
    if not any(v is not None for v in supplied) and not cache.manifest_path.exists():
        parser.error('authorized research inputs are required on first run; see docs/authorized_rnaseq.md')
    print(json.dumps(build(cache, *roots, cohort=cohort, clinical_file=args.clinical_file,
                           clinical_provenance=args.clinical_provenance), indent=2, sort_keys=True))
