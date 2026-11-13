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

from ds_bu_bioinformatics import cambridge
from ds_bu_bioinformatics import cambridge_clinical_source as publication
from ds_bu_bioinformatics.cambridge_pipeline import build, main
from ds_bu_bioinformatics.sources import SourceCache, sha256

FIXTURES = Path(__file__).parent / 'fixtures'


def publication_zip(conflict=False):
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    sheet = ET.Element('worksheet', xmlns=ns)
    data = ET.SubElement(sheet, 'sheetData')
    rows = [publication.HEADERS,
            ['TB01.0001_1_T', 'TUMOUR', 'RP', '10', 'Y', 'source BCR', '4'],
            ['TB01.0002_1_T', 'TUMOUR', 'RP', '12', 'N', 'N/A', 'N/A'],
            ['TB01.0003_1_T', 'TUMOUR', 'RP', 'N/A', 'N', 'N/A', 'N/A']]
    if conflict:
        rows.append(['TB01.0001_2_T', 'TUMOUR', 'RP', '10' if conflict == 'identical' else '20', 'Y', 'source BCR', '4'])
    for n, values in enumerate(rows, 1):
        row = ET.SubElement(data, 'row', r=str(n))
        for i, value in enumerate(values):
            cell = ET.SubElement(row, 'c', r=f'{chr(65+i)}{n}', t='inlineStr')
            ET.SubElement(ET.SubElement(cell, 'is'), 't').text = value
    xlsx = io.BytesIO()
    with zipfile.ZipFile(xlsx, 'w') as z:
        z.writestr('xl/worksheets/sheet1.xml', ET.tostring(sheet))
        z.writestr('xl/workbook.xml', f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Suppl. Table 2 clinical data" r:id="rId1"/></sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels', '<Relationships><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, 'w') as z:
        z.writestr('mmc2.xlsx', xlsx.getvalue())
    return archive.getvalue()


class CambridgeTest(unittest.TestCase):
    def setUp(self):
        self.metadata, self.expression, _ = cambridge.read_matrix(FIXTURES / 'cambridge_matrix.txt')
        self.clinical, self.samples, self.audit = cambridge.clinical_tables(self.metadata)

    def variant(self, old, new):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / 'matrix.txt'
            p.write_text((FIXTURES / 'cambridge_matrix.txt').read_text().replace(old, new))
            return cambridge.read_matrix(p)

    def test_selection_is_independent_of_missing_outcomes(self):
        self.assertEqual(self.samples.sample_id.tolist(), ['GSM1', 'GSM2', 'GSM3'])
        self.assertEqual(self.samples.patient_id.tolist(), ['TB01.0001', 'TB01.0002', 'TB01.0003'])
        self.assertEqual(self.audit.selection_reason.tolist(), ['primary_tumor'] * 3 + ['crpc', 'benign'])
        self.assertEqual(self.clinical.bcr_eligible.tolist(), [True, True, False])
        self.assertEqual(self.clinical.bcr_eligibility_reason.tolist(), ['eligible', 'eligible', 'missing_time'])
        self.assertTrue(pd.isna(self.clinical.bcr_months.iloc[2]))

    def test_event_and_censoring_use_different_time_fields(self):
        self.assertEqual(self.clinical.bcr_months.iloc[:2].tolist(), [4, 12])
        self.assertEqual(cambridge.decode_bcr('N', 'NA', '24'), (0, 24., 'eligible'))
        self.assertEqual(cambridge.decode_bcr('Y', None, '24'), (1, None, 'missing_time'))
        self.assertEqual(cambridge.decode_bcr('N/A', None, '24'), (None, None, 'missing_status'))
        self.assertEqual(cambridge.decode_bcr('Y', '0', '24'), (1, 0., 'zero_time_requires_review'))
        self.assertEqual(cambridge.decode_bcr('Y', '25', '24'),
                         (1, 25., 'event_after_recorded_followup_requires_review'))
        for args in [('X', '2', '24'), ('Y', '-1', '24'), ('Y', 'nan', '24'), ('N', None, 'inf')]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                cambridge.decode_bcr(*args)

    def test_preserve_features_drop_missing_across_all_primaries(self):
        selected = self.expression.loc[:, self.samples.sample_id]
        expression, genes, mapping = cambridge.map_features(FIXTURES / 'cambridge_platform.txt', selected)
        self.assertEqual(expression.index.tolist(), ['ILMN_1', 'ILMN_2', 'ILMN_4'])
        self.assertEqual(genes.gene_id.tolist(), ['10', '10', '30'])
        np.testing.assert_array_equal(expression.to_numpy(), [[2, 4, 6], [8, 10, 12], [5, 7, 9]])
        # Missing only in outcome-ineligible primary GSM3 still excludes the probe;
        # missing only in excluded CRPC GSM4 does not affect primary exports.
        self.assertEqual(mapping.loc[2, 'missing_selected_values'], 1)
        self.assertFalse(mapping.loc[2, 'selected'])
        self.assertTrue(mapping.loc[3, 'selected'])
        self.assertEqual(mapping.loc[5:6, 'mapping_status'].tolist(), ['gene_symbol_identity_conflict'] * 2)
        report = cambridge.validate(expression, self.clinical, self.samples, genes, mapping, self.audit)
        self.assertEqual(report['counts']['eligible_bcr_profiles'], 2)
        self.assertEqual(report['mapped_features_excluded_for_missing_expression'], 1)
        with self.assertRaisesRegex(ValueError, 'sample order'):
            cambridge.validate(expression, self.clinical, self.samples.iloc[::-1], genes, mapping, self.audit)

    def test_wrong_series_platform_truncation_and_invalid_values_rejected(self):
        for old, new in [('GSE70768', 'GSE70769'), ('GPL10558', 'GPL570'),
                         ('!series_matrix_table_end', ''), ('"ILMN_2"', '"ILMN_1"'),
                         ('"ILMN_1"\t2', '"ILMN_1"\tgarbage'), ('"ILMN_1"\t2', '"ILMN_1"\tinf')]:
            with self.subTest(old=old), self.assertRaises(ValueError):
                self.variant(old, new)

    def test_benign_title_variant_remains_excluded(self):
        meta = self.metadata.copy()
        meta.loc[4, 'title'] = 'benign tissue_robotic radical prostatectomy_TB01.0001'
        _, _, audit = cambridge.clinical_tables(meta)
        self.assertFalse(audit.selected.iloc[4])
        self.assertEqual(audit.selection_reason.iloc[4], 'benign')

    def test_unknown_specimen_and_duplicate_subject_rejected(self):
        bad = self.metadata.copy()
        bad.loc[0, 'title'] = 'unknown title'
        with self.assertRaisesRegex(ValueError, 'Unrecognized'):
            cambridge.clinical_tables(bad)
        bad = self.metadata.copy()
        bad.loc[1, 'title'] = bad.loc[0, 'title']
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            cambridge.clinical_tables(bad)

    def test_publication_joins_and_conflicting_subject_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'publication.zip'
            path.write_bytes(publication_zip())
            source = publication.read_clinical_source(path)
            # Publication provides follow-up missing from GEO, without overwriting GEO.
            source.loc[source.patient_id.eq('TB01.0003'), 'paper__followup_months'] = '10'
            clinical, _, _ = cambridge.clinical_tables(self.metadata, source)
            self.assertEqual(clinical.bcr_months.tolist(), [4, 12, 10])
            self.assertTrue(clinical.bcr_eligible.all())
            self.assertTrue(pd.isna(clinical.loc[2, 'geo__total_follow_up_months']))
            with self.assertRaisesRegex(ValueError, 'missing from publication'):
                cambridge.clinical_tables(self.metadata, source.iloc[:2])
            path.write_bytes(publication_zip(conflict=True))
            with self.assertRaisesRegex(ValueError, 'Conflicting publication'):
                publication.read_clinical_source(path)
            path.write_bytes(publication_zip(conflict='identical'))
            collapsed = publication.read_clinical_source(path)
            self.assertEqual(len(collapsed), 3)
            self.assertEqual(collapsed.paper__subject_id.iloc[0], 'TB01.0001_1_T;TB01.0001_2_T')

    def test_cli_rejects_unsafe_roots_and_snapshot(self):
        for args in [['--snapshot', '../other'], ['--output-root', 'data/raw/cambridge-v1'],
                     ['--output-root', 'data/anndata']]:
            with self.subTest(args=args), patch('urllib.request.urlopen', side_effect=AssertionError('CLI validation attempted network')), self.assertRaises(SystemExit):
                main(args)


class CambridgeExportTest(unittest.TestCase):
    def test_acquisition_csv_h5ad_and_offline_rebuild(self):
        class Response(io.BytesIO):
            headers = {}
            def __init__(self, content, url):
                super().__init__(content)
                self.url = url

        def transport(request, **kwargs):
            if request.full_url == cambridge.MATRIX_URL:
                return Response(gzip.compress((FIXTURES / 'cambridge_matrix.txt').read_bytes()), request.full_url)
            if request.full_url == cambridge.PLATFORM_URL:
                return Response((FIXTURES / 'cambridge_platform.txt').read_bytes(), request.full_url)
            if request.full_url == publication.URL:
                return Response(publication_zip(), request.full_url)
            raise AssertionError(request.full_url)

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cache = SourceCache(root / 'raw' / 'fixture', accessions=['GSE70768', 'GPL10558', 'PMC4588396'])
            roots = [root / name for name in ('processed', 'interim', 'anndata')]
            with patch('urllib.request.urlopen', side_effect=transport):
                report = build(cache, *roots)
            self.assertTrue(report['anndata_round_trip_valid'])
            self.assertTrue(report['export_round_trip_valid'])
            saved = ad.read_h5ad(roots[2] / 'fixture' / 'expression.h5ad')
            self.assertEqual(saved.shape, (3, 3))
            self.assertEqual(saved.var.gene_id.tolist(), ['10', '10', '30'])
            self.assertTrue(pd.isna(saved.obs.loc['GSM3', 'bcr_months']))
            self.assertFalse(saved.obs.loc['GSM3', 'bcr_eligible'])
            self.assertEqual(json.loads(saved.uns['metadata_json'])['geo_series'], 'GSE70768')
            rebuilt = [root / ('offline-' + name) for name in ('processed', 'interim', 'anndata')]
            with patch('urllib.request.urlopen', side_effect=AssertionError('unexpected network')):
                build(SourceCache(cache.root, offline=True, accessions=['GSE70768', 'GPL10558', 'PMC4588396']), *rebuilt)
            for original, copy in zip(roots, rebuilt):
                for p in (original / 'fixture').iterdir():
                    self.assertEqual(sha256(p), sha256(copy / 'fixture' / p.name), p.name)
            # A conflicting H5AD destination must prevent publishing either CSV or interim.
            blocked = [root / ('blocked-' + name) for name in ('processed', 'interim', 'anndata')]
            (blocked[2] / 'fixture').mkdir(parents=True)
            (blocked[2] / 'fixture' / 'keep.txt').write_text('existing')
            with self.assertRaisesRegex(ValueError, 'Existing outputs differ'):
                build(SourceCache(cache.root, offline=True, accessions=['GSE70768', 'GPL10558', 'PMC4588396']), *blocked)
            self.assertFalse((blocked[0] / 'fixture').exists())
            self.assertFalse((blocked[1] / 'fixture').exists())
            self.assertEqual((blocked[2] / 'fixture' / 'keep.txt').read_text(), 'existing')

    def test_wrong_cohort_cache_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            (p / 'manifest.json').write_text(json.dumps({'accessions': ['GSE21034'], 'sources': {}}))
            with self.assertRaisesRegex(ValueError, 'accessions differ'):
                SourceCache(p, accessions=['GSE70768', 'GPL10558', 'PMC4588396'])


if __name__ == '__main__':
    unittest.main()
