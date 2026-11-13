import json
import tempfile
import unittest
import urllib.parse
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.anndata_export import write_anndata
from ds_bu_bioinformatics.sources import SourceCache, sha256
from ds_bu_bioinformatics.tcga_pipeline import (
    DATAHUB_COMMIT, EXPR_URL, LEGACY_STUDY, STUDY, build, clean_expression, decode_endpoint, log2_expression,
    match_clinical, parse_number, read_export, read_rsem, select_samples, validate)

FIXTURES = Path(__file__).parent / 'fixtures'
CODES = {'0:CENSORED': 0, '1:PROGRESSION': 1}
SELECTED = ['TCGA-AA-0001-01', 'TCGA-AA-0002-01', 'TCGA-AA-0005-01']  # final cohort for the fixture


def clinical_records():
    """Synthetic cBioPortal records for the samples in fixtures/tcga_rsem.txt."""
    patients = {
        'TCGA-AA-0001': {'PFS_STATUS': '1:PROGRESSION', 'PFS_MONTHS': '12.5', 'OS_STATUS': '0:LIVING',
                         'OS_MONTHS': '20', 'AGE': '60', 'PATH_T_STAGE': 'T3A'},
        'TCGA-AA-0002': {'PFS_STATUS': '0:CENSORED', 'PFS_MONTHS': '30', 'AGE': '55', 'PATH_T_STAGE': 'T2C'},
        'TCGA-AA-0003': {'PFS_STATUS': '0:CENSORED', 'PFS_MONTHS': '10', 'AGE': '70'},
        'TCGA-AA-0005': {'PFS_STATUS': '0:CENSORED', 'PFS_MONTHS': '0', 'AGE': '[Not Available]'},
    }
    patient = [{'patientId': pid, 'clinicalAttributeId': attr, 'value': value, 'studyId': STUDY}
               for pid, attrs in patients.items() for attr, value in attrs.items()]
    sample_types = {'TCGA-AA-0001-01': 'Primary', 'TCGA-AA-0002-01': 'Primary',
                    'TCGA-AA-0003-11': 'Normal', 'TCGA-AA-0005-01': 'Primary'}
    sample = [{'sampleId': sid, 'patientId': sid[:12], 'clinicalAttributeId': 'SAMPLE_TYPE', 'value': value,
               'studyId': STUDY} for sid, value in sample_types.items()]
    # TCGA-AA-0004-01 is in the inventory but has no clinical records.
    inventory = [{'sampleId': sid, 'patientId': sid[:12], 'sampleType': 'Primary Solid Tumor', 'studyId': STUDY}
                 for sid in [*sample_types, 'TCGA-AA-0004-01']]
    legacy = [{'patientId': pid, 'clinicalAttributeId': attr, 'value': value, 'studyId': LEGACY_STUDY}
              for pid, attr, value in [('TCGA-AA-0001', 'GLEASON_SCORE', '7'), ('TCGA-AA-0002', 'GLEASON_SCORE', '9'),
                                       ('TCGA-AA-0001', 'GLEASON_PATTERN_PRIMARY', '4'),
                                       ('TCGA-AA-0001', 'OTHER_FIELD', 'ignored')]]
    return patient, sample, inventory, legacy


class FakeResponse:
    def __init__(self, payload, url):
        self.payload, self.url, self.headers = payload, url, {}

    def read(self, *args):
        payload, self.payload = self.payload, b''
        return payload

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def transport(request, **kwargs):
    url = request.full_url
    patient, sample, inventory, legacy = clinical_records()
    if url == EXPR_URL:
        return FakeResponse((FIXTURES / 'tcga_rsem.txt').read_bytes(), url)
    split = urllib.parse.urlsplit(url)
    query = urllib.parse.parse_qs(split.query)
    first = query.get('pageNumber') == ['0']
    if split.path.endswith(f'/studies/{STUDY}/clinical-attributes'):
        value = [{'clinicalAttributeId': k, 'patientAttribute': True, 'description': k}
                 for k in ['PFS_STATUS', 'PFS_MONTHS', 'OS_STATUS', 'OS_MONTHS']]
    elif split.path.endswith(f'/studies/{STUDY}/clinical-data'):
        value = (patient if query['clinicalDataType'] == ['PATIENT'] else sample) if first else []
    elif split.path.endswith(f'/studies/{STUDY}/samples'):
        value = inventory if first else []
    elif split.path.endswith(f'/studies/{LEGACY_STUDY}/clinical-data'):
        value = legacy if first else []
    else:
        raise AssertionError(url)
    return FakeResponse(json.dumps(value).encode(), url)


def tcga_cache(root, offline=False):
    return SourceCache(root, offline=offline, accessions=[STUDY, LEGACY_STUDY], user_agent='test')


class ExpressionTest(unittest.TestCase):
    def setUp(self):
        self.raw, self.features = read_rsem(FIXTURES / 'tcga_rsem.txt')

    def test_feature_ids_keep_repeated_entrez_rows_separate(self):
        self.assertEqual(self.features.feature_id.tolist(), ['354', '102', '103', '104_row4', '104_row5', '105', '106', '107'])
        self.assertEqual(self.features.gene_id.tolist()[3:5], ['104', '104'])
        self.assertEqual(self.features.gene_symbol.tolist()[3:5], ['OLDSYM', 'NEWSYM'])
        self.assertTrue(pd.isna(self.features.gene_symbol.iloc[5]))
        self.assertEqual(self.raw.shape, (8, 5))

    def test_missing_values_removed_negative_clipped_log2(self):
        expression, kept, audit, stats = clean_expression(self.raw, self.features, SELECTED)
        self.assertEqual(expression.columns.tolist(), SELECTED)
        self.assertNotIn('103', expression.index)
        self.assertEqual(audit.set_index('feature_id').loc['103', 'reason'], 'missing_values')
        self.assertEqual(stats['negative_rsem_values_clipped'], 1)
        self.assertEqual(kept.loc['106', 'TCGA-AA-0001-01'], -0.5)        # raw value preserved
        self.assertEqual(expression.loc['106', 'TCGA-AA-0001-01'], 0.0)   # log2(0 + 1)
        self.assertEqual(expression.loc['354', 'TCGA-AA-0001-01'], 10.0)  # log2(1023 + 1)
        np.testing.assert_array_equal(expression.loc['104_row4'], np.log2(8).round(4))

    def test_missing_value_filter_uses_selected_samples_only(self):
        # 107 is missing only in the normal sample TCGA-AA-0003-11, which is not selected.
        _, _, audit, stats = clean_expression(self.raw, self.features, SELECTED)
        audit = audit.set_index('feature_id')
        self.assertEqual(audit.loc['107', 'reason'], 'selected')
        self.assertEqual(audit.loc['103', 'missing_samples'], 1)
        self.assertEqual(stats['samples_with_missing_values'], 1)
        _, _, audit, _ = clean_expression(self.raw, self.features, ['TCGA-AA-0001-01', 'TCGA-AA-0005-01'])
        self.assertTrue(audit.selected.all())

    def read_variant(self, transform):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'rsem.txt'
            path.write_text(transform((FIXTURES / 'tcga_rsem.txt').read_text()))
            return read_rsem(path)

    def test_malformed_inputs_rejected(self):
        for name, transform in [('header', lambda s: s.replace('Entrez_Gene_Id', 'Entrez')),
                                ('barcode', lambda s: s.replace('TCGA-AA-0002-01', 'SAMPLE2')),
                                ('duplicate sample', lambda s: s.replace('TCGA-AA-0002-01', 'TCGA-AA-0001-01')),
                                ('entrez', lambda s: s.replace('GENEB\t102', 'GENEB\t')),
                                ('value', lambda s: s.replace('\t2047\t', '\tabc\t')),
                                ('infinite', lambda s: s.replace('\t2047\t', '\tinf\t'))]:
            with self.subTest(name), self.assertRaises(ValueError):
                self.read_variant(transform)


class SampleTest(unittest.TestCase):
    def test_barcode_rule_one_per_patient_and_qc_outlier(self):
        rng = np.random.default_rng(0)
        names = [f'TCGA-AA-{i:04d}-01' for i in range(14)] + ['TCGA-AA-0001-01B', 'TCGA-AA-0002-11']
        expression = pd.DataFrame(rng.normal(8, 0.1, (50, len(names))), columns=names)
        expression['TCGA-AA-0013-01'] = 0.0  # extreme on both median and genes detected
        audit = select_samples(expression, 3.0).set_index('sample_id').selection_reason
        self.assertEqual(audit['TCGA-AA-0002-11'], 'not_primary_tumor')
        self.assertEqual(audit['TCGA-AA-0001-01B'], 'extra_primary_for_patient')
        self.assertEqual(audit['TCGA-AA-0013-01'], 'qc_outlier')
        self.assertEqual((audit == 'primary_tumor').sum(), 13)

    def test_qc_ignores_features_missing_in_candidates(self):
        names = [f'TCGA-AA-{i:04d}-01' for i in range(10)]
        expression = pd.DataFrame(np.full((5, 10), 8.0), columns=names)
        expression.iloc[0, 0] = np.nan  # would change that sample's median if it were used
        audit = select_samples(expression, 3.0).set_index('sample_id')
        self.assertTrue(audit.median_log2.eq(8.0).all())
        self.assertTrue(audit.selected.all())


class ClinicalTest(unittest.TestCase):
    def setUp(self):
        raw, features = read_rsem(FIXTURES / 'tcga_rsem.txt')
        self.expression = log2_expression(raw)
        self.audit = select_samples(self.expression)
        self.records = clinical_records()

    def test_selection_prefixes_and_endpoints(self):
        clinical, samples, audit, _ = match_clinical(self.audit, *self.records)
        reasons = audit.set_index('sample_id').selection_reason.to_dict()
        self.assertEqual(reasons, {'TCGA-AA-0001-01': 'primary_tumor', 'TCGA-AA-0002-01': 'primary_tumor',
                                   'TCGA-AA-0003-11': 'not_primary_tumor', 'TCGA-AA-0004-01': 'no_clinical_record',
                                   'TCGA-AA-0005-01': 'primary_tumor'})
        self.assertEqual(clinical.sample_id.tolist(), samples.sample_id.tolist())
        row = clinical.set_index('sample_id')
        self.assertEqual(row.pfi_event.tolist(), [1, 0, 0])
        self.assertEqual(row.pfi_eligibility_reason.tolist(), ['eligible', 'eligible', 'zero_time_requires_review'])
        self.assertEqual(row.loc['TCGA-AA-0001-01', 'pfi_status'], 'PROGRESSION')
        self.assertEqual(row.os_eligibility_reason.tolist(), ['eligible', 'missing_status', 'missing_status'])
        self.assertEqual(row.loc['TCGA-AA-0001-01', 'path_t_group'], 'T3')
        self.assertEqual(row.loc['TCGA-AA-0002-01', 'gleason_score'], 9)
        self.assertTrue(pd.isna(row.loc['TCGA-AA-0005-01', 'age']))
        self.assertIn('legacy_patient__GLEASON_PATTERN_PRIMARY', clinical)
        self.assertNotIn('legacy_patient__OTHER_FIELD', clinical)
        self.assertIn('patient__PFS_STATUS', clinical)
        self.assertIn('sample__SAMPLE_TYPE', clinical)

    def test_invalid_age_or_gleason_fails(self):
        for attribute, value in [('AGE', 'sixty'), ('AGE', '-1'), ('GLEASON_SCORE', '7.5'), ('GLEASON_SCORE', '12')]:
            patient, sample, inventory, legacy = clinical_records()
            records = patient if attribute == 'AGE' else legacy
            next(r for r in records if r['clinicalAttributeId'] == attribute)['value'] = value
            with self.subTest(attribute=attribute, value=value), self.assertRaisesRegex(ValueError, 'Invalid'):
                match_clinical(self.audit, patient, sample, inventory, legacy)

    def test_parse_number_keeps_missing_tokens_missing(self):
        parsed = parse_number(pd.Series(['7', '[Not Available]', None, ' 8 ']), 'x', (2, 10), integer=True)
        self.assertEqual(parsed.tolist(), [7.0, pd.NA, pd.NA, 8.0])

    def test_inventory_disagreement_fails(self):
        patient, sample, inventory, legacy = self.records
        inventory[0]['patientId'] = 'TCGA-ZZ-9999'
        with self.assertRaisesRegex(ValueError, 'inventory'):
            match_clinical(self.audit, patient, sample, inventory, legacy)

    def test_non_primary_cbio_label_fails(self):
        patient, sample, inventory, legacy = self.records
        sample[0]['value'] = 'Metastatic'
        with self.assertRaisesRegex(ValueError, 'primary'):
            match_clinical(self.audit, patient, sample, inventory, legacy)

    def test_endpoint_decoding(self):
        self.assertEqual(decode_endpoint(None, '4', CODES), (None, 4.0, 'missing_status'))
        self.assertEqual(decode_endpoint('0:CENSORED', '[Not Available]', CODES), (0, None, 'missing_time'))
        self.assertEqual(decode_endpoint('0:CENSORED', '0', CODES), (0, 0.0, 'zero_time_requires_review'))
        self.assertEqual(decode_endpoint('1:PROGRESSION', '3.5', CODES), (1, 3.5, 'eligible'))
        for status, time in [('2:UNKNOWN', '2'), ('1:PROGRESSION', '-1'), ('1:PROGRESSION', 'nan')]:
            with self.subTest(status=status, time=time), self.assertRaises(ValueError):
                decode_endpoint(status, time, CODES)


class SourceTest(unittest.TestCase):
    def test_git_lfs_pointer_never_cached(self):
        pointer = b'version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 10\n'
        with tempfile.TemporaryDirectory() as temp:
            responses = [FakeResponse(pointer, 'https://example.test'), FakeResponse(b'real payload', 'https://example.test')]
            with patch('urllib.request.urlopen', side_effect=responses) as request, patch('time.sleep'):
                path = SourceCache(Path(temp)).fetch('data.txt', 'https://example.test')
            self.assertEqual(request.call_count, 2)
            self.assertEqual(path.read_bytes(), b'real payload')

    def test_anndata_layer_must_align(self):
        expression = pd.DataFrame([[1.0, 2.0]], index=pd.Index(['f1'], name='feature_id'), columns=['S1', 'S2'])
        clinical = pd.DataFrame({'sample_id': ['S1', 'S2']})
        genes = pd.DataFrame({'feature_id': ['f1'], 'gene_id': ['1']})
        with tempfile.TemporaryDirectory() as temp, self.assertRaisesRegex(ValueError, 'layer'):
            write_anndata(expression, clinical, genes, {}, Path(temp) / 'x.h5ad', layers={'rsem': expression[['S2', 'S1']]})


    def test_expression_url_pinned_to_datahub_commit(self):
        self.assertRegex(DATAHUB_COMMIT, r'^[0-9a-f]{40}$')
        self.assertIn(f'/datahub/{DATAHUB_COMMIT}/', EXPR_URL)


class EndToEndTest(unittest.TestCase):
    def test_download_export_offline_rebuild_and_anndata(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            modules = root / 'modules.csv'
            modules.write_text('gene,module\nKLK3,AR\nGENEC,AR\nMISSING,AR\nGENEB,Cycle\n')
            with patch('urllib.request.urlopen', side_effect=transport):
                report = build(tcga_cache(root / 'raw' / 'fixture'), root / 'out', z_threshold=3.0, module_file=modules)
            self.assertTrue(report['export_round_trip_valid'])
            self.assertTrue(report['anndata_round_trip_valid'])
            # TCGA-AA-0005-01 has PFI time 0 (ineligible), so it is not counted as censored.
            self.assertEqual(report['counts'], {'primary_profiles': 3, 'pfi_events': 1, 'pfi_censored': 1, 'gene_rows': 6})
            self.assertEqual(report['counts']['pfi_events'], report['endpoints']['pfi']['events'])
            self.assertEqual(report['eligible_pfi_profiles'], 2)
            self.assertEqual(report['retained_feature_rows'], 7)

            out = root / 'out' / 'fixture'
            expression, clinical, samples, genes, feature_audit, audit = read_export(out)
            self.assertEqual(expression.columns.tolist(), ['TCGA-AA-0001-01', 'TCGA-AA-0002-01', 'TCGA-AA-0005-01'])
            self.assertEqual(genes.feature_id.tolist(), ['354', '102', '104_row4', '104_row5', '105', '106', '107'])
            self.assertEqual(samples.sample_type_code.tolist(), ['01', '01', '01'])
            coverage = pd.read_csv(out / 'module_coverage.csv').set_index('module')
            self.assertEqual(coverage.loc['AR', ['genes', 'found', 'lost_to_missing_values']].tolist(), [3, 1, 1])
            for line in (out / 'checksums.sha256').read_text().splitlines():
                digest, name = line.split('  ', 1)
                self.assertEqual(sha256(out / name), digest, name)

            with patch('urllib.request.urlopen', side_effect=AssertionError('offline network attempted')):
                build(tcga_cache(root / 'raw' / 'fixture', offline=True), root / 'offline', root / 'offline-anndata',
                      module_file=modules)
            for path in out.iterdir():
                self.assertEqual(sha256(path), sha256(root / 'offline' / 'fixture' / path.name), path.name)

            import anndata as ad
            saved = ad.read_h5ad(root / 'anndata' / 'fixture' / 'expression.h5ad')
            self.assertEqual(saved.obs_names.tolist(), expression.columns.tolist())
            self.assertEqual(saved.var_names.tolist(), expression.index.tolist())
            np.testing.assert_array_equal(saved.X, expression.to_numpy().T)
            self.assertEqual(saved.layers['rsem'][0, saved.var_names.get_loc('106')], -0.5)
            self.assertEqual(json.loads(saved.uns['metadata_json'])['primary_endpoint'], 'pfi')

            with self.assertRaisesRegex(ValueError, 'order differs'):
                validate(expression, clinical, samples.iloc[::-1], genes, feature_audit, audit)

    def test_existing_different_outputs_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with patch('urllib.request.urlopen', side_effect=transport):
                build(tcga_cache(root / 'raw' / 'fixture'), root / 'out')
            with self.assertRaisesRegex(ValueError, 'differ'):
                build(tcga_cache(root / 'raw' / 'fixture', offline=True), root / 'out', z_threshold=1.0)


if __name__ == '__main__':
    unittest.main()
