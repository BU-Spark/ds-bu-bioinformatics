"""Acquire and export Cambridge GSE70768 primary tumors with BCR outcomes."""
from __future__ import annotations
from .provenance import code_hashes

import argparse
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from . import cambridge
from . import cambridge_clinical_source as clinical_source
from .anndata_export import write_anndata
from .expression import write_csv
from .sources import SourceCache, sha256, write_json
from .taylor_pipeline import check_destination, publish


def reload_csv(stage, expected):
    expression = pd.read_csv(stage / 'expression.csv.gz', index_col=0, dtype={'feature_id': str}, float_precision='round_trip')
    clinical = pd.read_csv(stage / 'clinical.csv', dtype='string', float_precision='round_trip')
    clinical['selected'] = clinical.selected.map({'True': True, 'False': False})
    clinical['bcr_event'] = clinical.bcr_event.astype('Int8')
    clinical['bcr_months'] = clinical.bcr_months.astype('Float64')
    clinical['bcr_eligible'] = clinical.bcr_eligible.map({'True': True, 'False': False})
    samples = pd.read_csv(stage / 'samples.csv', dtype='string')
    genes = pd.read_csv(stage / 'genes.csv', dtype='string')
    mapping = pd.read_csv(stage / 'feature_mapping.csv.gz', dtype={'feature_id': str, 'gene_id': str})
    audit = pd.read_csv(stage / 'sample_audit.csv', dtype={'sample_id': str})
    if expression.index.tolist() != expected.index.tolist() or expression.columns.tolist() != expected.columns.tolist() or not np.array_equal(expression.to_numpy(), expected.to_numpy()):
        raise ValueError('Cambridge CSV export changed expression')
    return cambridge.validate(expression, clinical, samples, genes, mapping, audit)


def build(cache, output_root: Path, interim_root: Path, anndata_root: Path, *, cohort=cambridge):
    cache.verify_all()
    metadata, source_expression, processing = cohort.read_matrix(cache.fetch(cohort.SERIES + '_series_matrix.txt.gz', cohort.MATRIX_URL))
    publication_path = cache.fetch('PMC4588396_supplementary.zip', clinical_source.URL)
    publication = clinical_source.read_clinical_source(publication_path, cohort=cohort.CLINICAL_COHORT)
    clinical, samples, audit = cohort.clinical_tables(metadata, publication)
    selected_source = source_expression.loc[:, samples.sample_id].copy()
    expression, genes, mapping = cambridge.map_features(cache.fetch('GPL10558.soft.txt', cambridge.PLATFORM_URL), selected_source)
    report = cambridge.validate(expression, clinical, samples, genes, mapping, audit)
    metadata_out = {
        'schema_version': 1, 'study': cohort.STUDY, 'geo_series': cohort.SERIES, 'platform': 'GPL10558',
        'snapshot': cache.root.name, 'source_manifest_sha256': sha256(cache.manifest_path),
        'sample_policy': cohort.SAMPLE_POLICY,
        'patient_id_provenance': 'Cohort subject label from GEO title, explicitly matched to one publication tumor subject row; no external registry linkage',
        'expression': {'orientation': 'features_by_samples', 'feature_key': 'Original ILMN probe ID',
            'scale': 'deposited processed signal; source reports log2 transform, quantile normalization and ComBat',
            'source_processing': processing, 'additional_transformations': 'none', 'aggregation': 'none',
            'missing_policy': 'Exclude entire probe if missing in ANY selected primary tumor, including tumors ineligible for BCR; record missing counts for every source feature; no imputation'},
        'mapping': {'source': 'Frozen GPL10558 Entrez_Gene_ID, Symbol, RefSeq_ID and Species',
            'policy': 'Require human species, one numeric Entrez ID and one symbol; quarantine inconsistent gene/symbol identities; keep separate ILMN features sharing the same gene; retain all feature decisions in audit',
            'mygene_refresh': False, 'symbol_validation': 'Deposited platform labels, not claimed current or HGNC-approved'},
        'bcr': {'endpoint': 'publication biochemical relapse (BCR)', 'source': 'PMC4588396 Supplementary Table 2, mmc2.xlsx',
            'source_url': clinical_source.URL, 'source_sha256': sha256(publication_path),
            'source_policy': 'Use publication status, event time and follow-up together for ALL selected primaries; preserve original GEO fields and publication subject IDs', 'status_codes': {'Y': 1, 'N': 0}, 'unit': 'months',
            'event_time_field': 'paper__bcr_months', 'censoring_time_field': 'paper__followup_months',
            'status_field': 'paper__bcr_status', 'missing_policy': 'Retain nullable event/time and explicit eligibility reason',
            'zero_time_policy': 'Retain but mark ineligible pending review',
            'time_conflict_policy': 'Retain BCR time and follow-up as deposited; flag event time later than recorded follow-up as ineligible pending source review',
            'definition_limit': 'No independently verified PSA threshold/confirmation or death handling; times interpreted as deposited postoperative BCR/follow-up durations'},
        'data_dictionary': {'expression.csv.gz': 'feature_id then sample_id columns; deposited values',
            'genes.csv': 'One row per expression feature; original ILMN feature_id, gene_id, gene_symbol, refseq_original',
            'clinical.csv': 'One row per primary tumor; source geo__ characteristics plus derived nullable bcr_* fields',
            'samples.csv': 'Expression sample order, GEO accession, title-derived subject ID and platform',
            'sample_audit.csv': 'Every GEO sample with title/type, selection and reason',
            'feature_mapping.csv.gz': 'Every expression feature in source order, annotation, mapping status, missing count, expression status and selected flag',
            'expression.h5ad': 'Separate AnnData root; samples x features, clinical in obs, genes in var, metadata_json in uns'},
        'limitations': ['Deposited annotation is not freshly revalidated with MyGene; differs from Taylor mapping source',
            'Missing-expression probes excluded without imputation; decisions are cohort-specific',
            'Subject identities are GEO-title-derived; no independent patient registry linkage',
            'Precise BCR criteria and death handling remain unverified; no raw assay QC'],
        'code_sha256': code_hashes(__name__, cohort)}
    stages = []
    try:
        for root in (output_root, interim_root, anndata_root):
            root.mkdir(parents=True, exist_ok=True)
            stages.append(Path(tempfile.mkdtemp(prefix='.' + cohort.CLINICAL_COHORT + '-', dir=root)))
        stage, interim_stage, ad_stage = stages
        write_csv(expression, stage / 'expression.csv.gz', index=True)
        for frame, name in [(clinical, 'clinical.csv'), (samples, 'samples.csv'), (genes, 'genes.csv'),
                            (mapping, 'feature_mapping.csv.gz'), (audit, 'sample_audit.csv')]:
            write_csv(frame, stage / name)
        write_csv(selected_source, interim_stage / 'feature_expression.csv.gz', index=True)
        reloaded = reload_csv(stage, expression)
        if reloaded['counts'] != report['counts']:
            raise ValueError('Cambridge CSV reload changed counts')
        write_anndata(expression, clinical, genes, metadata_out, ad_stage / 'expression.h5ad')
        report.update(export_round_trip_valid=True, anndata_round_trip_valid=True)
        write_json(stage / 'validation.json', report)
        for destination in (stage, ad_stage):
            write_json(destination / 'metadata.json', metadata_out)
            (destination / 'checksums.sha256').write_text(''.join(f'{sha256(p)}  {p.name}\n' for p in sorted(destination.iterdir())))
        for staged, root in zip(stages, (output_root, interim_root, anndata_root)):
            check_destination(staged, root / cache.root.name)
        for staged, root in zip(stages, (output_root, interim_root, anndata_root)):
            publish(staged, root / cache.root.name)
    finally:
        for staged in stages:
            shutil.rmtree(staged, ignore_errors=True)
    return report


def main(argv=None, *, cohort=cambridge):
    parser = argparse.ArgumentParser(description=f'Acquire and export {cohort.STUDY} {cohort.SERIES} primary tumors with BCR outcomes.')
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--snapshot', default=cohort.SNAPSHOT, help='Immutable cohort source snapshot name')
    parser.add_argument('--offline', action='store_true')
    for name in ('output-root', 'interim-root', 'anndata-root'):
        parser.add_argument('--' + name, type=Path)
    args = parser.parse_args(argv)
    if not args.snapshot or args.snapshot in {'.', '..'} or Path(args.snapshot).name != args.snapshot:
        parser.error('snapshot must be a single directory name')
    raw = args.data_dir / 'raw' / args.snapshot
    roots = [args.output_root or args.data_dir / 'processed', args.interim_root or args.data_dir / 'interim',
             args.anndata_root or args.data_dir / 'anndata']
    paths = [p.resolve() for p in [args.data_dir / 'raw', *roots]]
    if any(a.is_relative_to(b) or b.is_relative_to(a) for i, a in enumerate(paths) for b in paths[i+1:]):
        parser.error('raw, output, interim and AnnData directories must not overlap')
    cache = SourceCache(raw, offline=args.offline, accessions=[cohort.SERIES, 'GPL10558', 'PMC4588396'])
    print(json.dumps(build(cache, *roots, cohort=cohort), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
