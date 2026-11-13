import gzip
import io
import json
import tempfile
import unittest
import urllib.parse
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from ds_bu_bioinformatics.sources import SourceCache, GEO_URL, GPL_URL, sha256
from ds_bu_bioinformatics.expression import read_geo
from ds_bu_bioinformatics.annotation import parse_platform, resolve_mappings, select_features
from ds_bu_bioinformatics.clinical import match_clinical, decode_dfs, pivot_clinical, paged
from ds_bu_bioinformatics.taylor_pipeline import build, publish
from ds_bu_bioinformatics.validation import validate
from ds_bu_bioinformatics.anndata_export import write_anndata

FIXTURES = Path(__file__).parent / 'fixtures'


def hit(query, gene, symbol):
    return {'query': query, 'entrezgene': gene, 'symbol': symbol, 'taxid': 9606, '_score': 10, 'refseq': {'rna': [query + '.1']}}


def clinical_records():
    patient = []
    for pid, status, months in [('PAT1', '1:Recurred', '12'), ('PAT3', '0:DiseaseFree', '24')]:
        for attr, value in [('DFS_STATUS', status), ('DFS_MONTHS', months), ('GLEASON_SCORE', '7')]:
            patient.append({'patientId': pid, 'clinicalAttributeId': attr, 'value': value, 'studyId': 'prad_mskcc'})
    sample = [{'sampleId': sid, 'patientId': pid, 'clinicalAttributeId': 'SAMPLE_TYPE', 'value': 'Primary', 'studyId': 'prad_mskcc'}
              for sid, pid in [('PCA3', 'PAT3'), ('PCA1', 'PAT1')]]
    inventory = [{'sampleId': r['sampleId'], 'patientId': r['patientId']} for r in sample]
    return patient, sample, inventory


class ExpressionTest(unittest.TestCase):
    def read_variant(self, transform):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / 'matrix.txt'
            p.write_text(transform((FIXTURES / 'matrix.txt').read_text()))
            return read_geo(p)

    def test_source_and_empty_characteristic(self):
        meta, exp, processing = read_geo(FIXTURES / 'matrix.txt')
        self.assertEqual(exp.shape, (4, 3))
        self.assertTrue(pd.isna(meta.loc[1, 'tumor_type']))
        self.assertEqual(exp.loc['1', 'GSM1'], 2)
        self.assertIn('quantile normalization', processing[0])

    def test_truncation_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Truncated|Malformed'):
            self.read_variant(lambda s: s.replace('!series_matrix_table_end', ''))

    def test_duplicate_feature_rejected(self):
        with self.assertRaisesRegex(ValueError, 'duplicate feature'):
            self.read_variant(lambda s: s.replace('"2"\t6', '"1"\t6'))

    def test_missing_expression_rejected(self):
        with self.assertRaisesRegex(ValueError, 'non-finite'):
            self.read_variant(lambda s: s.replace('"2"\t6', '"2"\tnan'))

    def test_wrong_series_rejected(self):
        with self.assertRaisesRegex(ValueError, 'series'):
            self.read_variant(lambda s: s.replace('GSE21034', 'GSE99999'))

    def test_wrong_platform_rejected(self):
        with self.assertRaisesRegex(ValueError, 'platform'):
            self.read_variant(lambda s: s.replace('GPL10264', 'GPL999'))

    def test_column_order_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError, 'columns'):
            self.read_variant(lambda s: s.replace('"ID_REF"\t"GSM1"\t"GSM2"', '"ID_REF"\t"GSM2"\t"GSM1"'))


class AnnotationTest(unittest.TestCase):
    def setUp(self):
        self.platform = parse_platform(FIXTURES / 'platform.txt')
        self.payload = [hit('NM_1', 10, 'A'), hit('NM_2', 10, 'A'), hit('NM_3', 20, 'B'),
                        hit('NM_3', 30, 'C'), {'query': 'NM_4', 'notfound': True}]

    def test_ambiguity_excluded_and_duplicate_genes_kept_separate(self):
        mapping = resolve_mappings(self.platform, self.payload, ['1', '2', '3', '4'])
        self.assertEqual(mapping.loc[mapping.feature_id.eq('3'), 'mapping_status'].tolist(), ['ambiguous', 'ambiguous'])
        self.assertEqual(mapping.loc[mapping.feature_id.eq('1'), 'refseq_original'].iloc[0], 'NM_1.2')
        _, exp, _ = read_geo(FIXTURES / 'matrix.txt')
        gene, genes = select_features(exp[['GSM1', 'GSM3']], mapping)
        np.testing.assert_array_equal(gene.values, [[2, 8], [6, 4]])
        self.assertEqual(genes.feature_id.tolist(), ['1', '2'])
        self.assertEqual(genes.gene_id.tolist(), ['10', '10'])
        self.assertEqual(genes.gene_symbol.tolist(), ['A', 'A'])
        self.assertEqual(mapping.feature_id.nunique(), 4)

    def test_annotation_loss_rejected(self):
        with self.assertRaisesRegex(ValueError, 'sets differ'):
            resolve_mappings(self.platform, self.payload, ['1'])

    def test_symbol_collision_quarantined(self):
        payload = [hit('NM_1', 10, 'A'), hit('NM_2', 20, 'A')]
        result = resolve_mappings(self.platform, payload, ['1', '2', '3', '4'])
        self.assertFalse(result.selected.any())
        self.assertIn('gene_symbol_identity_conflict', result.mapping_status.tolist())

    def test_numeric_string_gene_id_and_exact_refseq(self):
        payload = [hit('NM_1', '10', 'A'), dict(hit('NM_2', '20', 'B'), refseq={'rna': ['NM_999.1']})]
        result = resolve_mappings(self.platform, payload, ['1', '2', '3', '4'])
        self.assertTrue(result.loc[result.feature_id.eq('1'), 'selected'].iloc[0])
        self.assertEqual(result.loc[result.feature_id.eq('2'), 'mapping_status'].iloc[0], 'refseq_mismatch')

    def test_nonhuman_invalid_candidate(self):
        payload = [dict(hit('NM_1', 10, 'A'), taxid=10090)]
        result = resolve_mappings(self.platform, payload, ['1', '2', '3', '4'])
        self.assertFalse(result.selected.any())


class ClinicalTest(unittest.TestCase):
    def setUp(self):
        self.meta, _, _ = read_geo(FIXTURES / 'matrix.txt')
        self.patient, self.sample, self.inventory = clinical_records()

    def test_shuffled_records_match_real_patient_ids(self):
        clinical, audit, _ = match_clinical(self.meta, self.patient[::-1], self.sample, self.inventory)
        self.assertEqual(clinical.sample_id.tolist(), ['PCA1', 'PCA3'])
        self.assertEqual(clinical.patient_id.tolist(), ['PAT1', 'PAT3'])
        self.assertEqual(clinical.dfs_event.tolist(), [1, 0])
        self.assertEqual(audit.selected.tolist(), [True, False, True])

    def test_missing_join_fails(self):
        with self.assertRaisesRegex(ValueError, 'Unmatched'):
            match_clinical(self.meta, self.patient, self.sample[:1], self.inventory)

    def test_patient_inventory_disagreement_fails(self):
        self.inventory[0]['patientId'] = 'BAD'
        with self.assertRaisesRegex(ValueError, 'inventory'):
            match_clinical(self.meta, self.patient, self.sample, self.inventory)

    def test_multiple_primaries_per_patient_fails(self):
        self.sample[0]['patientId'] = 'PAT1'
        self.inventory[0]['patientId'] = 'PAT1'
        with self.assertRaisesRegex(ValueError, 'Multiple primary'):
            match_clinical(self.meta, self.patient, self.sample, self.inventory)

    def test_conflicting_duplicate_fails(self):
        records = [*self.patient, dict(self.patient[0], value='0:DiseaseFree')]
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            pivot_clinical(records, 'PATIENT')

    def test_identical_duplicate_counted(self):
        _, count = pivot_clinical([*self.patient, self.patient[0]], 'PATIENT')
        self.assertEqual(count, 1)

    def test_missing_and_zero_outcomes_preserved(self):
        self.assertEqual(decode_dfs(None, '4'), (None, 4.0, 'missing_status'))
        self.assertEqual(decode_dfs('0:DiseaseFree', None), (0, None, 'missing_time'))
        self.assertEqual(decode_dfs('0:DiseaseFree', '0'), (0, 0.0, 'zero_time_requires_review'))

    def test_invalid_outcomes_fail(self):
        for status, time in [('2:Unknown', '2'), ('1:Recurred', '-1'), ('1:Recurred', 'nan')]:
            with self.subTest(status=status, time=time), self.assertRaises(ValueError):
                decode_dfs(status, time)


class FakeResponse(io.BytesIO):
    def __init__(self, payload, url, headers=None):
        super().__init__(payload)
        self.url = url
        self.headers = headers or {}


class SourceTest(unittest.TestCase):
    def test_empty_response_retried_automatically(self):
        with tempfile.TemporaryDirectory() as temp:
            responses = [FakeResponse(b'', 'https://example.test'), FakeResponse(b'complete', 'https://example.test')]
            with patch('urllib.request.urlopen', side_effect=responses) as request, patch('time.sleep'):
                path = SourceCache(Path(temp)).fetch('input', 'https://example.test')
            self.assertEqual(request.call_count, 2)
            self.assertEqual(path.read_bytes(), b'complete')

    def test_content_length_mismatch_retried(self):
        with tempfile.TemporaryDirectory() as temp:
            responses = [FakeResponse(b'abc', 'https://example.test', {'Content-Length': '20'}),
                         FakeResponse(b'complete', 'https://example.test', {'Content-Length': '8'})]
            with patch('urllib.request.urlopen', side_effect=responses) as request, patch('time.sleep'):
                path = SourceCache(Path(temp)).fetch('input', 'https://example.test')
            self.assertEqual(request.call_count, 2)
            self.assertEqual(path.read_bytes(), b'complete')

    def test_html_error_page_never_cached_as_json(self):
        with tempfile.TemporaryDirectory() as temp:
            responses = [FakeResponse(b'<html>Temporary failure</html>', 'https://example.test'),
                         FakeResponse(b'[]', 'https://example.test')]
            with patch('urllib.request.urlopen', side_effect=responses), patch('time.sleep'):
                path = SourceCache(Path(temp)).fetch('input.json', 'https://example.test')
            self.assertEqual(json.loads(path.read_text()), [])

    def test_truncated_gzip_retried_before_caching(self):
        with tempfile.TemporaryDirectory() as temp:
            valid = gzip.compress(b'complete input')
            responses = [FakeResponse(valid[:-8], 'https://example.test'), FakeResponse(valid, 'https://example.test')]
            with patch('urllib.request.urlopen', side_effect=responses) as request, patch('time.sleep'):
                path = SourceCache(Path(temp)).fetch('input.gz', 'https://example.test')
            self.assertEqual(request.call_count, 2)
            self.assertEqual(gzip.decompress(path.read_bytes()), b'complete input')

    def test_persistent_empty_response_leaves_no_cache(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch('urllib.request.urlopen', side_effect=lambda *a, **k: FakeResponse(b'', 'https://example.test')) as request, patch('time.sleep'):
                with self.assertRaisesRegex(ValueError, 'Empty HTTP response'):
                    SourceCache(Path(temp)).fetch('input', 'https://example.test')
            self.assertEqual(request.call_count, 3)
            self.assertEqual(list(Path(temp).iterdir()), [])

    def test_cache_reuse_corruption_and_offline_absence(self):
        with tempfile.TemporaryDirectory() as temp:
            cache = SourceCache(Path(temp))
            with patch('urllib.request.urlopen', return_value=FakeResponse(b'abc', 'https://example.test')) as request:
                path = cache.fetch('input', 'https://example.test')
                cache.fetch('input', 'https://example.test')
                self.assertEqual(request.call_count, 1)
            with patch('urllib.request.urlopen', side_effect=AssertionError('network')):
                offline = SourceCache(Path(temp), offline=True)
                offline.fetch('input', 'https://example.test')
                with self.assertRaisesRegex(ValueError, 'Offline'):
                    offline.fetch('absent', 'https://example.test')
                path.write_bytes(b'corrupt')
                with self.assertRaisesRegex(ValueError, 'corrupt'):
                    offline.fetch('input', 'https://example.test')

    def test_partial_download_not_published(self):
        class BrokenResponse(FakeResponse):
            def read(self, *args):
                raise ConnectionError('interrupted')
        with tempfile.TemporaryDirectory() as temp:
            with patch('urllib.request.urlopen', side_effect=lambda *a, **k: BrokenResponse(b'x', 'https://example.test')), patch('time.sleep'):
                with self.assertRaises(ConnectionError):
                    SourceCache(Path(temp)).fetch('input', 'https://example.test')
            self.assertFalse((Path(temp) / 'input').exists())
            self.assertFalse((Path(temp) / 'input.part').exists())

    def test_repeated_pagination_rejected(self):
        class Cache:
            manifest = {'sources': {f'p_{i:04d}.json': {'response_metadata': {}} for i in range(2)}}
            def json(self, *args):
                return [{'sampleId': 'A'}]
        with self.assertRaisesRegex(ValueError, 'Repeated'):
            paged(Cache(), 'samples', 'p', page_size=1)

    def test_incomplete_pagination_rejected(self):
        class Cache:
            manifest = {'sources': {f'p_{i:04d}.json': {'response_metadata': {'total-count': '2'}} for i in range(2)}}
            def json(self, name, *args):
                return [{'sampleId': 'A'}] if '0000' in name else []
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            paged(Cache(), 'samples', 'p', page_size=1)


class EndToEndTest(unittest.TestCase):
    def test_download_export_and_offline_byte_rebuild(self):
        patient, sample, inventory = clinical_records()
        def transport(request, **kwargs):
            url = request.full_url
            if url == GEO_URL:
                return FakeResponse(gzip.compress((FIXTURES / 'matrix.txt').read_bytes(), mtime=0), url)
            if url == GPL_URL:
                return FakeResponse((FIXTURES / 'platform.txt').read_bytes(), url)
            if url.endswith('/clinical-attributes'):
                value = [{'clinicalAttributeId': k, 'patientAttribute': True, 'description': k} for k in ['DFS_STATUS', 'DFS_MONTHS']]
            elif '/clinical-data?' in url:
                q = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
                value = (patient if q['clinicalDataType'] == ['PATIENT'] else sample) if q['pageNumber'] == ['0'] else []
            elif '/samples?' in url:
                q = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
                value = inventory if q['pageNumber'] == ['0'] else []
            elif url.endswith('/metadata'):
                value = {'build_version': 'test'}
            elif url.endswith('/query'):
                value = [hit('NM_1', 10, 'A'), hit('NM_2', 10, 'A'), hit('NM_3', 20, 'B'), {'query': 'NM_4', 'notfound': True}]
            else:
                raise AssertionError(url)
            return FakeResponse(json.dumps(value).encode(), url)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with patch('urllib.request.urlopen', side_effect=transport):
                report = build(SourceCache(root / 'raw' / 'fixture'), root / 'out', root / 'interim')
            self.assertTrue(report['export_round_trip_valid'])
            self.assertTrue(report['anndata_round_trip_valid'])
            self.assertEqual(report['counts']['primary_profiles'], 2)  # Never force 131.
            with patch('urllib.request.urlopen', side_effect=AssertionError('offline network attempted')):
                build(SourceCache(root / 'raw' / 'fixture', offline=True), root / 'offline', root / 'offline-interim', root / 'offline-anndata')
            metadata_path = root / 'out' / 'fixture' / 'metadata.json'
            original_metadata = metadata_path.read_bytes()
            offline_cache = SourceCache(root / 'raw' / 'fixture', offline=True)
            with patch('ds_bu_bioinformatics.taylor_pipeline.sha256',
                       side_effect=lambda path: '0' * 64 if path.suffix == '.py' else sha256(path)):
                build(offline_cache, root / 'out', root / 'interim')
                self.assertEqual(metadata_path.read_bytes(), original_metadata)
                samples_path = metadata_path.with_name('samples.csv')
                original_samples = samples_path.read_bytes()
                samples_path.write_text('corrupt')
                with self.assertRaisesRegex(ValueError, 'Existing outputs differ'):
                    build(offline_cache, root / 'out', root / 'interim')
                self.assertEqual(samples_path.read_text(), 'corrupt')
                samples_path.write_bytes(original_samples)
            for path in (root / 'out' / 'fixture').iterdir():
                self.assertEqual(sha256(path), sha256(root / 'offline' / 'fixture' / path.name), path.name)
            genes = pd.read_csv(root / 'out' / 'fixture' / 'expression.csv.gz', index_col=0)
            np.testing.assert_array_equal(genes.values, [[2, 8], [6, 4], [9, 3]])
            import anndata as ad
            saved = ad.read_h5ad(root / 'anndata' / 'fixture' / 'expression.h5ad')
            self.assertEqual(saved.obs_names.tolist(), ['PCA1', 'PCA3'])
            self.assertEqual(saved.var_names.tolist(), ['1', '2', '3'])
            self.assertEqual(saved.var.gene_symbol.tolist(), ['A', 'A', 'B'])
            np.testing.assert_array_equal(saved.X, genes.to_numpy().T)
            rebuilt = ad.read_h5ad(root / 'offline-anndata' / 'fixture' / 'expression.h5ad')
            np.testing.assert_array_equal(rebuilt.X, saved.X)
            pd.testing.assert_frame_equal(rebuilt.obs, saved.obs)
            pd.testing.assert_frame_equal(rebuilt.var, saved.var)
            with self.assertRaisesRegex(ValueError, 'order differs'):
                clinical = pd.read_csv(root / 'out' / 'fixture' / 'clinical.csv')
                samples = pd.read_csv(root / 'out' / 'fixture' / 'samples.csv').iloc[::-1]
                validate(genes, clinical, samples, pd.DataFrame(), pd.DataFrame(), pd.DataFrame())

    def test_h5ad_preserves_missing_clinical_values(self):
        import anndata as ad
        expression = pd.DataFrame([[2., 6.], [8., 4.]],
                                  index=pd.Index(['f2', 'f1'], name='feature_id'), columns=['S2', 'S1'])
        clinical = pd.DataFrame({'sample_id': ['S1', 'S2'],
                                 'label': pd.Series([None, 'primary'], dtype='string'),
                                 'dfs_event': pd.Series([None, 1], dtype='Int64'),
                                 'dfs_months': pd.Series([None, 12.5], dtype='Float64'),
                                 'dfs_eligible': pd.Series([False, True], dtype='boolean')})
        genes = pd.DataFrame({'feature_id': ['f1', 'f2'], 'gene_id': ['10', '10'], 'gene_symbol': ['A', 'A']})
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'expression.h5ad'
            write_anndata(expression, clinical, genes, {'aggregation': 'none'}, path)
            saved = ad.read_h5ad(path)
            self.assertEqual(saved.obs_names.tolist(), ['S2', 'S1'])
            self.assertEqual(saved.var_names.tolist(), ['f2', 'f1'])
            self.assertTrue(saved.obs.loc['S1', ['label', 'dfs_event', 'dfs_months']].isna().all())
            self.assertEqual(saved.obs.loc['S2', 'dfs_months'], 12.5)
            self.assertEqual(saved.var.gene_symbol.tolist(), ['A', 'A'])

    def test_existing_different_outputs_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); stage = root / 'stage'; dest = root / 'dest'
            stage.mkdir(); dest.mkdir(); (stage / 'x').write_text('new'); (dest / 'x').write_text('old')
            with self.assertRaisesRegex(ValueError, 'differ'):
                publish(stage, dest)
            self.assertEqual((dest / 'x').read_text(), 'old')


if __name__ == '__main__':
    unittest.main()
