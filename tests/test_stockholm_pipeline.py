import gzip
import io
import json
import tempfile
import unittest
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

from ds_bu_bioinformatics import cambridge, stockholm
from ds_bu_bioinformatics import cambridge_clinical_source as publication
from ds_bu_bioinformatics.sources import SourceCache, sha256
from ds_bu_bioinformatics.stockholm_pipeline import build, main
from test_cambridge_pipeline import publication_zip

FIXTURES = Path(__file__).parent / 'fixtures'


def stockholm_publication_zip():
    # Add Stockholm to the synthetic two-cohort publication without changing
    # Cambridge records; this tests that each adapter selects its own subjects.
    with zipfile.ZipFile(io.BytesIO(publication_zip())) as archive:
        with zipfile.ZipFile(io.BytesIO(archive.read('mmc2.xlsx'))) as book:
            members = {name: book.read(name) for name in book.namelist()}
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    sheet = ET.fromstring(members['xl/worksheets/sheet1.xml'])
    data = sheet.find(f'{{{ns}}}sheetData')
    records = [
        ['STKHLM1', 'TUMOUR', 'RP', '10', 'Y ', 'FORMAL (>0.2)', '4'],
        ['STKHLM2', 'TUMOUR', 'RP', '12', 'N', 'N/A', 'N/A'],
        ['STKHLM3', 'TUMOUR', 'RP', 'UNKNOWN', 'N', 'N/A', 'UNKNOWN'],
        ['STKHLM4', 'TUMOUR', 'RP', 'N/A', 'N/A', 'N/A', 'N/A']]
    for n, values in enumerate(records, 5):
        row = ET.SubElement(data, f'{{{ns}}}row', r=str(n))
        for i, value in enumerate(values):
            cell = ET.SubElement(row, f'{{{ns}}}c', r=f'{chr(65+i)}{n}', t='inlineStr')
            ET.SubElement(ET.SubElement(cell, f'{{{ns}}}is'), f'{{{ns}}}t').text = value
    members['xl/worksheets/sheet1.xml'] = ET.tostring(sheet)
    xlsx = io.BytesIO()
    with zipfile.ZipFile(xlsx, 'w') as book:
        for name, content in members.items():
            book.writestr(name, content)
    result = io.BytesIO()
    with zipfile.ZipFile(result, 'w') as archive:
        archive.writestr('mmc2.xlsx', xlsx.getvalue())
    return result.getvalue()


class StockholmTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        path = Path(self.temp.name) / 'publication.zip'
        path.write_bytes(stockholm_publication_zip())
        self.publication = publication.read_clinical_source(path, cohort='stockholm')
        self.metadata, self.expression, _ = stockholm.read_matrix(FIXTURES / 'stockholm_matrix.txt')
        self.clinical, self.samples, self.audit = stockholm.clinical_tables(self.metadata, self.publication)

    def test_cohort_specific_subject_selection_and_missing_outcomes(self):
        self.assertEqual(self.publication.patient_id.tolist(), ['STKHLM1', 'STKHLM2', 'STKHLM3', 'STKHLM4'])
        self.assertEqual(self.samples.patient_id.tolist(), self.publication.patient_id.tolist())
        self.assertTrue(self.audit.selected.all())
        self.assertEqual(self.clinical.bcr_eligible.tolist(), [True, True, False, False])
        self.assertEqual(self.clinical.bcr_eligibility_reason.tolist(), ['eligible', 'eligible', 'missing_time', 'missing_status'])
        self.assertEqual(self.clinical.bcr_months.iloc[:2].tolist(), [4, 12])
        self.assertTrue(pd.isna(self.clinical.bcr_event.iloc[3]))
        # Cambridge remains isolated when reading the same workbook.
        path = Path(self.temp.name) / 'publication.zip'
        cambridge_subjects = publication.read_clinical_source(path)
        self.assertTrue(cambridge_subjects.patient_id.str.startswith('TB').all())

    def test_finite_probe_policy_keeps_all_samples_and_separate_gene_rows(self):
        exp, genes, mapping = cambridge.map_features(FIXTURES / 'cambridge_platform.txt', self.expression)
        self.assertEqual(exp.shape, (3, 4))
        self.assertEqual(genes.gene_id.tolist(), ['10', '10', '30'])
        np.testing.assert_array_equal(exp.iloc[:2].to_numpy(), [[2, 4, 6, 8], [8, 10, 12, 14]])
        report = cambridge.validate(exp, self.clinical, self.samples, genes, mapping, self.audit)
        self.assertEqual(report['counts']['primary_profiles'], 4)
        self.assertEqual(report['counts']['eligible_bcr_profiles'], 2)
        self.assertEqual(report['mapped_features_excluded_for_missing_expression'], 1)
        # Missing value in an ineligible subject still excludes the whole probe.
        self.assertFalse(mapping.loc[mapping.feature_id.eq('ILMN_3'), 'selected'].iloc[0])

    def test_specimen_identity_and_publication_coverage_fail_closed(self):
        for title in ['matched benign tissue_robotic radical prostatectomy_STKHLM1',
                      'tumour tissue_robotic radical prostatectomy_TB01.0001']:
            bad = self.metadata.copy()
            bad.loc[0, 'title'] = title
            with self.subTest(title=title), self.assertRaises(ValueError):
                stockholm.clinical_tables(bad, self.publication)
        bad = self.metadata.copy()
        bad.loc[1, 'title'] = bad.loc[0, 'title']
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            stockholm.clinical_tables(bad, self.publication)
        with self.assertRaisesRegex(ValueError, 'missing from publication'):
            stockholm.clinical_tables(self.metadata, self.publication.iloc[:3])
        bad = self.publication.copy()
        bad.loc[0, 'paper__specimen_type'] = 'Biopsy'
        with self.assertRaisesRegex(ValueError, 'primary RP'):
            stockholm.clinical_tables(self.metadata, bad)

    def test_reader_rejects_cambridge_and_cli_validates_roots_without_network(self):
        with self.assertRaisesRegex(ValueError, 'GSE70769'):
            stockholm.read_matrix(FIXTURES / 'cambridge_matrix.txt')
        for args in [['--snapshot', '../other'], ['--output-root', 'data/raw/stockholm-v1'],
                     ['--output-root', 'data/anndata']]:
            with self.subTest(args=args), patch('urllib.request.urlopen', side_effect=AssertionError('network')), self.assertRaises(SystemExit):
                main(args)

    def test_online_acquisition_csv_h5ad_and_identical_offline_rebuild(self):
        class Response(io.BytesIO):
            headers = {}
            def __init__(self, value, url):
                super().__init__(value)
                self.url = url

        def transport(request, **kwargs):
            url = request.full_url
            if url == stockholm.MATRIX_URL:
                payload = gzip.compress((FIXTURES / 'stockholm_matrix.txt').read_bytes())
            elif url == cambridge.PLATFORM_URL:
                payload = (FIXTURES / 'cambridge_platform.txt').read_bytes()
            elif url == publication.URL:
                payload = stockholm_publication_zip()
            else:
                raise AssertionError(url)
            return Response(payload, url)

        root = Path(self.temp.name)
        accessions = ['GSE70769', 'GPL10558', 'PMC4588396']
        cache = SourceCache(root / 'raw' / 'fixture', accessions=accessions)
        roots = [root / name for name in ('processed', 'interim', 'anndata')]
        with patch('urllib.request.urlopen', side_effect=transport):
            report = build(cache, *roots)
        self.assertTrue(report['export_round_trip_valid'])
        self.assertTrue(report['anndata_round_trip_valid'])
        data = ad.read_h5ad(roots[2] / 'fixture' / 'expression.h5ad')
        self.assertEqual(data.shape, (4, 3))
        self.assertTrue(pd.isna(data.obs.loc['GSM4', 'bcr_event']))
        self.assertTrue(pd.isna(data.obs.loc['GSM3', 'bcr_months']))
        self.assertEqual(data.var.gene_id.tolist(), ['10', '10', '30'])
        metadata = json.loads(data.uns['metadata_json'])
        self.assertEqual(metadata['geo_series'], 'GSE70769')
        self.assertEqual(metadata['study'], 'Stockholm')
        self.assertIn('STKHLM', metadata['sample_policy'])
        copies = [root / ('offline-' + name) for name in ('processed', 'interim', 'anndata')]
        with patch('urllib.request.urlopen', side_effect=AssertionError('offline network')):
            build(SourceCache(cache.root, offline=True, accessions=accessions), *copies)
        for original, copied in zip(roots, copies):
            for p in (original / 'fixture').iterdir():
                self.assertEqual(sha256(p), sha256(copied / 'fixture' / p.name), p.name)


if __name__ == '__main__':
    unittest.main()
