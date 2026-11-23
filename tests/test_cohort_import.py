import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from ds_bu_bioinformatics.cohort_import import BundleCohort, COHORTS, main
from ds_bu_bioinformatics.geo_cohort_pipeline import build
from ds_bu_bioinformatics.output_check import check_dataset
from ds_bu_bioinformatics.sources import SourceCache, sha256


def bundle(root, cohort='erho'):
    root.mkdir()
    data = {
        'expression': 'feature_id,GSM1,GSM2,GSM3\np1,2,3,4\np2,4,6,8\np3,1,,3\n',
        'genes': 'feature_id,gene_id,gene_symbol\np1,10,A\np2,10,A\np3,20,B\n',
        'samples': 'sample_id,patient_id,platform_id\nGSM1,P1,GPL5188\nGSM2,P2,GPL5188\n',
        'clinical': 'sample_id,metastasis_event,metastasis_time,metastasis_unit\nGSM2,0,2,years\nGSM1,1,,months\n',
        'sample_audit': 'sample_id,selected,selection_reason\nGSM1,True,primary\nGSM2,True,primary\nGSM3,False,benign\n'}
    spec = {'schema_version': 1, 'study_accession': COHORTS[cohort][0], 'organism': 'Homo sapiens',
        'platforms': ['GPL5188'], 'expression_kind': 'log_microarray', 'expression_scale': 'Unchanged supplied normalized signal',
        'selection_policy': 'Fixture primary samples, exclude benign; exact source patient mapping',
        'annotation_source': 'Synthetic platform fixture', 'source_description': 'Synthetic research-input fixture',
        'primary_endpoint': 'metastasis', 'endpoints': {'metastasis': {'definition': 'Source metastatic event',
            'time_origin': 'surgery', 'censoring': 'last metastasis-free follow-up', 'sample_linkage': 'GSM column to fixture record'}},
        'files': {}}
    for logical, content in data.items():
        path = root / (logical + ('.csv.gz' if logical == 'expression' else '.csv'))
        if logical == 'expression':
            path.write_bytes(gzip.compress(content.encode(), mtime=0))
        else:
            path.write_text(content)
        spec['files'][logical] = {'path': path.name, 'sha256': sha256(path), 'source': f'fixture:{logical}'}
    manifest = root / 'manifest.json'
    manifest.write_text(json.dumps(spec))
    return manifest, spec


class CohortImportTest(unittest.TestCase):
    def test_all_cohort_identities_have_a_strict_local_import_path(self):
        for key in COHORTS:
            with self.subTest(cohort=key), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                manifest, _ = bundle(root / 'supplied', key)
                cohort = BundleCohort(key, manifest)
                cache = SourceCache(root / 'raw' / cohort.SNAPSHOT, accessions=cohort.ACCESSIONS)
                with patch('urllib.request.urlopen', side_effect=AssertionError('unexpected HTTP')):
                    c, s, a, x, g, *_ = cohort.prepare(cache)
                self.assertEqual(c.sample_id.tolist(), ['GSM1', 'GSM2'])
                self.assertEqual(c.metastasis_event.tolist(), [1, 0])
                self.assertEqual(c.metastasis_eligible.tolist(), [False, True])
                self.assertEqual(c.metastasis_months.iloc[1], 24)
                self.assertEqual(x.columns.tolist(), s.sample_id.tolist())
                self.assertEqual(g.gene_id.tolist(), ['10', '10', '20'])

    def test_import_passes_output_checker_and_replays_offline_without_input_bundle(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, _ = bundle(root / 'supplied')
            cohort = BundleCohort('erho', manifest)
            cache = SourceCache(root / 'raw' / cohort.SNAPSHOT, accessions=cohort.ACCESSIONS)
            roots = [root / p for p in ('processed', 'interim', 'anndata')]
            report = build(cache, *roots, cohort=cohort)
            self.assertEqual(report['counts']['retained_feature_rows'], 2)
            self.assertEqual(check_dataset(roots[0] / cohort.SNAPSHOT, roots[2] / cohort.SNAPSHOT)['errors'], [])
            # Removing the original directory is safe: every input is frozen in raw/.
            for path in manifest.parent.iterdir():
                path.unlink()
            manifest.parent.rmdir()
            copies = [root / ('offline-' + p) for p in ('processed', 'interim', 'anndata')]
            with patch('urllib.request.urlopen', side_effect=AssertionError('offline HTTP')):
                replay = BundleCohort('erho')
                build(SourceCache(cache.root, offline=True, accessions=replay.ACCESSIONS), *copies, cohort=replay)
            for source, copy in zip(roots, copies):
                for path in (source / cohort.SNAPSHOT).iterdir():
                    self.assertEqual(sha256(path), sha256(copy / cohort.SNAPSHOT / path.name), path.name)

    def test_wrong_cohort_checksum_path_or_endpoint_provenance_fail(self):
        for field, value in [('study_accession', 'GSE_OTHER'), ('organism', 'Mus musculus'),
                             ('primary_endpoint', 'os'), ('selection_policy', '')]:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                manifest, spec = bundle(root / 'supplied')
                spec[field] = value
                manifest.write_text(json.dumps(spec))
                cohort = BundleCohort('erho', manifest)
                with self.assertRaises(ValueError):
                    cohort.prepare(SourceCache(root / 'raw' / 'test', accessions=cohort.ACCESSIONS))
        for key, value in [('sha256', '0' * 64), ('path', '../outside.csv')]:
            with self.subTest(key=key), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                manifest, spec = bundle(root / 'supplied')
                spec['files']['clinical'][key] = value
                manifest.write_text(json.dumps(spec))
                cohort = BundleCohort('erho', manifest)
                with self.assertRaises(ValueError):
                    cohort.prepare(SourceCache(root / 'raw' / 'test', accessions=cohort.ACCESSIONS))

    def test_duplicate_patient_and_incomplete_audit_are_rejected(self):
        for logical, old, new in [('samples', 'GSM2,P2', 'GSM2,P1'),
                                  ('sample_audit', 'GSM3,False,benign\n', ''),
                                  ('clinical', 'GSM2,0,2,years', 'GSM999,0,2,years')]:
            with self.subTest(logical=logical), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                manifest, spec = bundle(root / 'supplied')
                path = manifest.parent / spec['files'][logical]['path']
                path.write_text(path.read_text().replace(old, new))
                spec['files'][logical]['sha256'] = sha256(path)
                manifest.write_text(json.dumps(spec))
                cohort = BundleCohort('erho', manifest)
                with self.assertRaises(ValueError):
                    cohort.prepare(SourceCache(root / 'raw' / 'test', accessions=cohort.ACCESSIONS))

    def test_raw_counts_freeze_library_totals_before_annotation_filtering(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, spec = bundle(root / 'supplied', 'moscow')
            spec['expression_kind'] = 'raw_counts'
            path = manifest.parent / spec['files']['expression']['path']
            path.write_bytes(gzip.compress(b'feature_id,GSM1,GSM2,GSM3\np1,2,3,4\np2,4,6,8\np3,100,200,300\n', mtime=0))
            spec['files']['expression']['sha256'] = sha256(path)
            genes = manifest.parent / spec['files']['genes']['path']
            genes.write_text(genes.read_text().replace('p3,20,B', 'p3,,---'))
            spec['files']['genes']['sha256'] = sha256(genes)
            manifest.write_text(json.dumps(spec))
            cohort = BundleCohort('moscow', manifest)
            cache = SourceCache(root / 'raw' / cohort.SNAPSHOT, accessions=cohort.ACCESSIONS)
            build(cache, root / 'processed', root / 'interim', root / 'anndata', cohort=cohort)
            from ds_bu_bioinformatics.module_scores import load_dataset
            expression, *_, transform = load_dataset(root / 'processed' / cohort.SNAPSHOT)
            import numpy as np
            np.testing.assert_allclose(expression.loc['p1'], np.log2(np.array([2 / 106, 3 / 209]) * 1e6 + 1))
            self.assertIn('log2(CPM + 1)', transform)
            metadata = json.loads((root / 'processed' / cohort.SNAPSHOT / 'metadata.json').read_text())
            self.assertEqual(metadata['expression']['library_sizes'], {'GSM1': 106, 'GSM2': 209})
            raw = pd.read_csv(root / 'processed' / cohort.SNAPSHOT / 'expression.csv.gz', index_col=0)
            self.assertEqual(raw.loc['p1'].tolist(), [2, 3])

    def test_raw_count_missingness_and_unknown_scale_fail(self):
        for kind in ('raw_counts', 'unknown'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                manifest, spec = bundle(root / 'supplied')
                spec['expression_kind'] = kind
                manifest.write_text(json.dumps(spec))
                cohort = BundleCohort('erho', manifest)
                with self.assertRaises(ValueError):
                    cohort.prepare(SourceCache(root / 'raw' / 'test', accessions=cohort.ACCESSIONS))

    def test_duplicate_expression_columns_are_not_silently_renamed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, spec = bundle(root / 'supplied')
            path = manifest.parent / spec['files']['expression']['path']
            path.write_bytes(gzip.compress(b'feature_id,GSM1,GSM1\np1,1,2\n', mtime=0))
            spec['files']['expression']['sha256'] = sha256(path)
            manifest.write_text(json.dumps(spec))
            cohort = BundleCohort('erho', manifest)
            with self.assertRaisesRegex(ValueError, 'unique'):
                cohort.prepare(SourceCache(root / 'raw' / 'test', accessions=cohort.ACCESSIONS))


if __name__ == '__main__':
    unittest.main()
