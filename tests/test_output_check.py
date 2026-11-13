import json
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.anndata_export import write_anndata
from ds_bu_bioinformatics.expression import write_csv
from ds_bu_bioinformatics.output_check import check_dataset, main
from ds_bu_bioinformatics.sources import sha256, write_json


def write_checksums(folder):
    (folder / 'checksums.sha256').write_text(''.join(
        f'{sha256(p)}  {p.name}\n' for p in sorted(folder.iterdir()) if p.name != 'checksums.sha256'))


def make_dataset(data_dir, name='demo-v1', scale=1.0, edit=None):
    """Write a small dataset in the shared output format. edit(tables) may change tables before writing."""
    samples_ids = ['S1', 'S2', 'S3']
    expression = pd.DataFrame(np.array([[5.0, 6.0, 7.0], [1.5, 2.5, 3.5]]) * scale,
                              index=pd.Index(['f1', 'f2'], name='feature_id'), columns=samples_ids)
    tables = {
        'samples.csv': pd.DataFrame({'sample_id': samples_ids, 'patient_id': ['P1', 'P2', 'P3'], 'platform_id': 'GPLX'}),
        'clinical.csv': pd.DataFrame({'sample_id': samples_ids, 'bcr_event': pd.array([1, 0, None], dtype='Int8'),
                                      'bcr_months': pd.array([10.0, 20.0, None], dtype='Float64'),
                                      'bcr_eligible': [True, True, False],
                                      'bcr_eligibility_reason': ['eligible', 'eligible', 'missing_status'],
                                      'bcr_time_unit': 'months'}),
        'genes.csv': pd.DataFrame({'feature_id': ['f1', 'f2'], 'gene_id': ['1', '2'], 'gene_symbol': ['A', 'B']}),
        'sample_audit.csv': pd.DataFrame({'sample_id': [*samples_ids, 'S4'], 'selected': [True, True, True, False],
                                          'selection_reason': ['primary_tumor'] * 3 + ['not_primary_tumor']}),
        'feature_audit.csv': pd.DataFrame({'feature_id': ['f1', 'f2'], 'selected': True}),
    }
    metadata = {'snapshot': name}
    validation = {'structurally_valid': True, 'export_round_trip_valid': True, 'anndata_round_trip_valid': True}
    if edit:
        edit(tables)
    processed, anndata_dir = data_dir / 'processed' / name, data_dir / 'anndata' / name
    processed.mkdir(parents=True)
    anndata_dir.mkdir(parents=True)
    write_csv(expression, processed / 'expression.csv.gz', index=True)
    for file_name, frame in tables.items():
        write_csv(frame, processed / file_name)
    write_json(processed / 'metadata.json', metadata)
    write_json(processed / 'validation.json', validation)
    write_anndata(expression, tables['clinical.csv'], tables['genes.csv'], metadata, anndata_dir / 'expression.h5ad')
    write_json(anndata_dir / 'metadata.json', metadata)
    write_checksums(processed)
    write_checksums(anndata_dir)
    return processed, anndata_dir


class OutputCheckTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.data = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def check(self, **kwargs):
        return check_dataset(*make_dataset(self.data, **kwargs))

    def test_valid_dataset_passes_with_counts_from_files(self):
        result = self.check()
        self.assertEqual(result['errors'], [])
        self.assertEqual(result['warnings'], [])
        self.assertEqual((result['samples'], result['feature_rows'], result['unique_genes']), (3, 2, 2))
        self.assertEqual(result['endpoints'], {'bcr': {'eligible': 2, 'events': 1, 'time_unit': 'months'}})

    def test_order_and_identity_problems_fail(self):
        def reorder_clinical(tables):
            tables['clinical.csv'] = tables['clinical.csv'].iloc[::-1]

        def repeat_patient(tables):
            tables['samples.csv'].loc[2, 'patient_id'] = 'P1'

        def reorder_genes(tables):
            tables['genes.csv'] = tables['genes.csv'].iloc[::-1]

        for edit, message in [(reorder_clinical, 'clinical.csv: sample order'),
                              (repeat_patient, 'one sample per patient'),
                              (reorder_genes, 'genes.csv: feature order')]:
            with self.subTest(message), tempfile.TemporaryDirectory() as temp:
                errors = check_dataset(*make_dataset(Path(temp), edit=edit))['errors']
                self.assertTrue(any(message in e for e in errors), errors)

    def test_endpoint_problems_fail(self):
        def eligible_without_time(tables):
            tables['clinical.csv']['bcr_months'] = pd.array([10.0, None, None], dtype='Float64')

        def missing_unit(tables):
            tables['clinical.csv'] = tables['clinical.csv'].drop(columns='bcr_time_unit')

        def bad_event(tables):
            tables['clinical.csv']['bcr_event'] = pd.array([1, 2, None], dtype='Int8')

        for edit, message in [(eligible_without_time, 'positive time'), (missing_unit, 'lacks'),
                              (bad_event, 'other than 0, 1')]:
            with self.subTest(message), tempfile.TemporaryDirectory() as temp:
                errors = check_dataset(*make_dataset(Path(temp), edit=edit))['errors']
                self.assertTrue(any(message in e for e in errors), errors)

    def test_checksum_and_missing_file_fail(self):
        processed, anndata_dir = make_dataset(self.data)
        (processed / 'clinical.csv').write_text((processed / 'clinical.csv').read_text().replace('months', 'Months', 1))
        self.assertTrue(any('checksum mismatch' in e for e in check_dataset(processed, anndata_dir)['errors']))
        (anndata_dir / 'expression.h5ad').unlink()
        self.assertTrue(any('missing expression.h5ad' in e for e in check_dataset(processed, anndata_dir)['errors']))

    def test_anndata_must_match_csv(self):
        processed, anndata_dir = make_dataset(self.data)
        other = self.data / 'other'
        _, other_anndata = make_dataset(other, scale=2.0)
        (anndata_dir / 'expression.h5ad').write_bytes((other_anndata / 'expression.h5ad').read_bytes())
        write_checksums(anndata_dir)
        errors = check_dataset(processed, anndata_dir)['errors']
        self.assertTrue(any('X differs' in e for e in errors), errors)

    def test_unlogged_expression_warns_but_passes(self):
        result = self.check(scale=1000.0)
        self.assertEqual(result['errors'], [])
        self.assertTrue(any('log-transformed' in w for w in result['warnings']))

    def test_cli_checks_every_snapshot_and_sets_exit_code(self):
        make_dataset(self.data, 'good-v1')
        output = StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(['--data-dir', str(self.data), '--json', str(self.data / 'report.json')]), 0)
        self.assertIn('| good-v1 | PASS | 3 |', output.getvalue())
        self.assertEqual(json.loads((self.data / 'report.json').read_text())[0]['dataset'], 'good-v1')
        make_dataset(self.data, 'bad-v1', edit=lambda t: t.update({'genes.csv': t['genes.csv'].iloc[::-1]}))
        with redirect_stdout(StringIO()) as output:
            self.assertEqual(main(['--data-dir', str(self.data)]), 1)
        self.assertIn('| bad-v1 | FAIL', output.getvalue())


if __name__ == '__main__':
    unittest.main()
