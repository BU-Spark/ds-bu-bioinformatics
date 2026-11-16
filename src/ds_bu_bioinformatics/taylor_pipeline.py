"""Download, audit, and export Taylor/MSKCC GSE21034 GPL10264 + prad_mskcc."""
from __future__ import annotations
import argparse
import json
import shutil
import tempfile
from pathlib import Path
import pandas as pd
from .sources import SourceCache, GEO_URL, GPL_URL, CBIO_BASE, sha256, write_json
from .expression import read_geo, write_csv
from .annotation import parse_platform, fetch_mappings, resolve_mappings, select_features
from .clinical import paged, match_clinical, DFS_CODES
from .validation import validate, validate_export
from .anndata_export import write_anndata


def check_destination(stage: Path, destination: Path):
    if destination.exists():
        old = {str(p.relative_to(destination)): sha256(p) for p in destination.rglob('*') if p.is_file()}
        new = {str(p.relative_to(stage)): sha256(p) for p in stage.rglob('*') if p.is_file()}
        if old != new:
            raise ValueError(f'Existing outputs differ: {destination}; choose a new output root or snapshot')


def publish(stage: Path, destination: Path):
    check_destination(stage, destination)
    if destination.exists():
        shutil.rmtree(stage)
    else:
        stage.replace(destination)


def build(cache, output_root: Path, interim_root: Path, anndata_root: Path | None = None):
    anndata_root = anndata_root if anndata_root is not None else output_root.parent / 'anndata'
    cache.verify_all()
    matrix_path = cache.fetch('GSE21034-GPL10264_series_matrix.txt.gz', GEO_URL)
    metadata, source_expression, processing = read_geo(matrix_path)
    platform = parse_platform(cache.fetch('GPL10264.soft.txt', GPL_URL))
    attributes = cache.json('cbio_clinical_attributes.json', f'{CBIO_BASE}/clinical-attributes')
    patient = paged(cache, 'clinical-data', 'cbio_patient_clinical', {'clinicalDataType': 'PATIENT', 'projection': 'DETAILED'})
    sample = paged(cache, 'clinical-data', 'cbio_sample_clinical', {'clinicalDataType': 'SAMPLE', 'projection': 'DETAILED'})
    inventory = paged(cache, 'samples', 'cbio_samples', {'projection': 'DETAILED'})
    clinical, audit, duplicate_stats = match_clinical(metadata, patient, sample, inventory)
    # Use source matrix order explicitly; never rely on merge or API response order.
    order = metadata.loc[metadata.tumor_type.eq('Primary tumor'), 'sample_id'].tolist()
    clinical = clinical.set_index('sample_id').loc[order].reset_index()
    samples = clinical[['sample_id', 'geo_accession', 'cbio_sample_id', 'patient_id', 'platform_id', 'tumor_type']].copy()
    features = source_expression.loc[:, samples.geo_accession].copy()
    features.columns = pd.Index(samples.sample_id, name=None)
    cache.json('mygene/metadata.json', 'https://mygene.info/v3/metadata')
    payload = fetch_mappings(cache, platform)
    mapping = resolve_mappings(platform, payload, features.index)
    expression, genes = select_features(features, mapping)
    checks = validate(expression, clinical, samples, genes, mapping, audit)
    checks.update(duplicate_stats)
    attribute_map = {a['clinicalAttributeId']: a for a in attributes if a.get('patientAttribute')}
    if not {'DFS_STATUS', 'DFS_MONTHS'}.issubset(attribute_map):
        raise ValueError('DFS attribute definitions missing')
    # Evidence is raw GEO processing text and downloaded study-specific attribute definitions.
    metadata_out = {
        'schema_version': 2, 'study': 'Taylor/MSKCC', 'geo_series': 'GSE21034', 'platform': 'GPL10264', 'cbio_study': 'prad_mskcc',
        'snapshot': cache.root.name, 'source_manifest_sha256': sha256(cache.manifest_path),
        'expression': {'orientation': 'features_by_samples', 'feature_key': 'original GPL10264 feature ID; gene IDs and symbols may repeat',
            'source_processing': processing, 'scale': 'deposited_signal; exact log base not established',
            'normalization': 'GEO reports RMA background adjustment and quantile normalization',
            'additional_transformations': 'none', 'aggregation': 'none; retain each unambiguously mapped feature separately',
            'scale_evidence': 'GSE21034-GPL10264_series_matrix.txt.gz !Sample_data_processing',
            'log2_assumption_permitted': False},
        'mapping': {'source': 'MyGene.info v3; frozen complete batch responses',
            'policy': 'Human candidates with exact version-stripped RefSeq RNA match and one consistent Entrez/symbol identity; ambiguity excluded; no highest-score selection',
            'refseq_policy': 'Preserve original accession; strip terminal numeric version only for lookup',
            'symbol_validation': 'MyGene labels; not claimed HGNC-approved'},
        'dfs': {'endpoint': 'DFS', 'status_codes': DFS_CODES, 'unit': 'months',
            'definitions': {k: attribute_map[k] for k in ['DFS_STATUS', 'DFS_MONTHS']},
            'time_origin': 'initial treatment per cBio study attributes; specific treatment/surgery origin unverified',
            'recurrence_criteria': 'precise PSA threshold, confirmation rule and death handling unverified',
            'zero_time_policy': 'retain but exclude from survival eligibility pending source review',
            'missing_policy': 'retain nullable event/time and eligibility reason; never assume censored'},
        'clinical_attributes': attributes,
        'clinical_field_provenance': {'patient__*': 'PATIENT clinical records, pivoted by actual patientId',
            'sample__*': 'SAMPLE clinical records, pivoted by sampleId',
            'dfs_event': 'explicit code map from patient__DFS_STATUS',
            'dfs_months': 'numeric patient__DFS_MONTHS; no unit conversion'},
        'sample_policy': 'GEO Primary tumor; require cBio Primary label when supplied; exact sample ID and inventory patient association; reject multiple primaries per patient',
        'data_dictionary': {
            'expression.csv.gz': 'feature_id followed by sample_id columns; deposited numeric signal',
            'samples.csv': 'One row per expression column; GEO and actual cBio sample/patient identifiers',
            'clinical.csv': 'Same sample order; original prefixed covariates plus nullable DFS/event/time and eligibility',
            'genes.csv': 'One row per expression feature; feature_id, gene_id, gene_symbol and RefSeq accessions',
            'feature_mapping.csv.gz': 'Every source feature and every returned candidate; original/query RefSeq, score, mapping status and selection reason',
            'sample_audit.csv': 'Every GEO profile, selection reason, available cBio identifiers',
            'expression.h5ad': 'Separate AnnData folder: samples × features; clinical in obs, gene/RefSeq annotations in var, metadata_json in uns'},
        'limitations': ['Exact deposited log base remains unverified; no extra transformation applied',
            'DFS is not relabeled BCR; precise recurrence definition remains unverified',
            'No raw CEL quality control; this workflow starts with deposited processed expression'],
        'code_sha256': {p.name: sha256(p) for p in sorted(Path(__file__).parent.glob('*.py'))}}
    metadata_out['output_schemas'] = {name: {column: {'dtype': str(frame[column].dtype), 'nullable': bool(frame[column].isna().any())} for column in frame.columns} for name, frame in [('clinical.csv', clinical), ('samples.csv', samples), ('genes.csv', genes), ('feature_mapping.csv.gz', mapping), ('sample_audit.csv', audit)]}
    # Rebuild with the published provenance when verifying the existing pair.
    # Every staged file still has to match its published counterpart byte for byte.
    published_metadata = output_root / cache.root.name / 'metadata.json'
    if published_metadata.is_file() and (anndata_root / cache.root.name).is_dir():
        try:
            published = json.loads(published_metadata.read_text())
        except ValueError:
            published = None  # The destination check will report the mismatch.
        if isinstance(published, dict) and isinstance(published.get('code_sha256'), dict):
            metadata_out['code_sha256'] = published['code_sha256']
    checks['expression_scale_verified'] = False
    checks['precise_dfs_definition_verified'] = False
    # Deterministic outputs allow byte-for-byte offline verification, including gzip files.
    output_root.mkdir(parents=True, exist_ok=True)
    interim_root.mkdir(parents=True, exist_ok=True)
    anndata_root.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.taylor-', dir=output_root))
    interim_stage = Path(tempfile.mkdtemp(prefix='.taylor-', dir=interim_root))
    anndata_stage = Path(tempfile.mkdtemp(prefix='.taylor-', dir=anndata_root))
    try:
        write_csv(expression, stage / 'expression.csv.gz', index=True)
        for frame, name in [(clinical, 'clinical.csv'), (samples, 'samples.csv'), (genes, 'genes.csv'),
                            (mapping, 'feature_mapping.csv.gz'), (audit, 'sample_audit.csv')]:
            write_csv(frame, stage / name)
        write_csv(features, interim_stage / 'feature_expression.csv.gz', index=True)
        reloaded = validate_export(stage, expression)
        if reloaded['counts'] != checks['counts']:
            raise ValueError('Export round-trip changed counts')
        checks['export_round_trip_valid'] = True
        write_anndata(expression, clinical, genes, metadata_out, anndata_stage / 'expression.h5ad')
        checks['anndata_round_trip_valid'] = True
        write_json(anndata_stage / 'metadata.json', metadata_out)
        (anndata_stage / 'checksums.sha256').write_text(''.join(
            f'{sha256(p)}  {p.name}\n' for p in sorted(anndata_stage.iterdir())))
        write_json(stage / 'metadata.json', metadata_out)
        write_json(stage / 'validation.json', checks)
        checksums = '\n'.join(f'{sha256(p)}  {p.name}' for p in sorted(stage.iterdir())) + '\n'
        (stage / 'checksums.sha256').write_text(checksums)
        # Check every destination before publishing any output format.
        for staged, root in [(interim_stage, interim_root), (stage, output_root), (anndata_stage, anndata_root)]:
            check_destination(staged, root / cache.root.name)
        # Publish only after every export has validated.
        publish(interim_stage, interim_root / cache.root.name)
        publish(stage, output_root / cache.root.name)
        publish(anndata_stage, anndata_root / cache.root.name)
    finally:
        shutil.rmtree(stage, ignore_errors=True)
        shutil.rmtree(interim_stage, ignore_errors=True)
        shutil.rmtree(anndata_stage, ignore_errors=True)
    return checks


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--snapshot', default='taylor-v2', help='Immutable input snapshot name; use a new name to refresh sources')
    parser.add_argument('--offline', action='store_true', help='Forbid HTTP; require all inputs and matching checksums')
    parser.add_argument('--output-root', type=Path, help='Alternate processed root, useful for offline verification')
    parser.add_argument('--interim-root', type=Path, help='Alternate intermediate root')
    parser.add_argument('--anndata-root', type=Path, help='Alternate AnnData output root (default: data/anndata)')
    args = parser.parse_args(argv)
    if not args.snapshot or args.snapshot in {'.', '..'} or Path(args.snapshot).name != args.snapshot:
        parser.error('snapshot must be a single directory name')
    cache = SourceCache(args.data_dir / 'raw' / args.snapshot, offline=args.offline)
    output_root = args.output_root or args.data_dir / 'processed'
    interim_root = args.interim_root or args.data_dir / 'interim'
    anndata_root = args.anndata_root or args.data_dir / 'anndata'
    paths = [cache.root.resolve(), output_root.resolve(), interim_root.resolve(), anndata_root.resolve()]
    if any(a.is_relative_to(b) or b.is_relative_to(a) for i, a in enumerate(paths) for b in paths[i+1:]):
        parser.error('raw, output, interim and AnnData directories must not overlap')
    report = build(cache, output_root, interim_root, anndata_root)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
