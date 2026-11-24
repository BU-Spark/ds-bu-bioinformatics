import contextlib
import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from ds_bu_bioinformatics import nakagawa, nakagawa_pipeline
from ds_bu_bioinformatics.cambridge import read_matrix
from ds_bu_bioinformatics.geo_cohort_pipeline import finite_features
from ds_bu_bioinformatics.output_check import check_dataset
from ds_bu_bioinformatics.sources import MYGENE_URL, SourceCache, sha256
from test_mortensen_long import matrix_bytes


def inputs():
    payloads = {}
    for platform, start, order in [('GPL5858', 0, [1, 2, 3]), ('GPL5873', 3, [3, 1, 2])]:
        records = []
        for i, subject in enumerate(order):
            records.append({'sample_id': f'GSM{start + i + 1}',
                'title': f'Patient {subject:03d} Core Set (custom panel)' if platform == 'GPL5873' else f'Patient {subject:03d} Core Set (cancer panel)',
                'Case-Control Group': {1: 'SYS', 2: 'PSA', 3: 'NED'}[subject],
                'Prostate cancer specific death': 'Yes' if subject == 1 else 'No',
                'Years of Followup after RRP': str(subject + 4),
                'Time (yrs) from PSA recurrence to systemic disease': '10' if subject == 2 else 'NA'})
        rows = [['GI_1-S', *[2 * p for p in order]], ['GI_2-S', *[3 * p for p in order]],
                ['unknown', *[4 * p for p in order]]]
        payloads[nakagawa.MATRIX_URLS[platform]] = gzip.compress(matrix_bytes(nakagawa.SERIES, platform, records, rows), mtime=0)
        payloads[nakagawa.PLATFORM_URLS[platform]] = ('\n'.join([f'^PLATFORM = {platform}', '!platform_table_begin',
            'ID\tSequenceSource\tGB_ACC', 'GI_1-S\tRefSeq\tNM_1.1', 'GI_2-S\tRefSeq\tNM_2.2',
            'unknown\tGenBank\tAI1234', '!platform_table_end']) + '\n').encode()
    payloads[MYGENE_URL] = json.dumps([{'query': f'NM_{n}', 'symbol': symbol, 'entrezgene': n * 10,
        'taxid': 9606, 'refseq': {'rna': f'NM_{n}.1'}} for n, symbol in [(1, 'A'), (2, 'B')]]).encode()
    return payloads


def transport(payloads):
    class Response(io.BytesIO):
        headers = {}
        def __init__(self, value, url):
            super().__init__(value)
            self.url = url
    return patch('urllib.request.urlopen', side_effect=lambda request, **kwargs: Response(payloads[request.full_url], request.full_url))


class NakagawaTest(unittest.TestCase):
    def prepare(self, root):
        cache = SourceCache(root / 'raw' / nakagawa.SNAPSHOT, accessions=nakagawa.ACCESSIONS)
        with transport(inputs()), contextlib.redirect_stdout(io.StringIO()):
            prepared = nakagawa.prepare(cache)
        return cache, prepared

    def test_exact_two_panel_join_reorders_assays_and_preserves_panel_features(self):
        with tempfile.TemporaryDirectory() as temp:
            _, (clinical, samples, audit, source, annotations, extra, *_rest) = self.prepare(Path(temp))
        self.assertEqual(samples.patient_id.tolist(), ['001', '002', '003'])
        self.assertEqual(samples.gpl5873_assay_id.tolist(), ['GSM5', 'GSM6', 'GSM4'])
        self.assertEqual(source.index.tolist(), [p + ':' + f for p in nakagawa.PANELS for f in ['GI_1-S', 'GI_2-S', 'unknown']])
        np.testing.assert_array_equal(source.iloc[[0, 3]].to_numpy(), [[2, 4, 6], [2, 4, 6]])
        expression, genes, mapping = finite_features(source, annotations)
        self.assertEqual(genes.gene_id.tolist(), ['10', '20', '10', '20'])
        self.assertEqual(expression.shape, (4, 3))
        self.assertEqual(len(extra['assay_linkage.csv']), 6)
        self.assertEqual(audit.sample_id.tolist(), samples.sample_id.tolist())

    def test_group_labels_do_not_become_systemic_events_and_followup_is_separate(self):
        with tempfile.TemporaryDirectory() as temp:
            _, (clinical, *_rest) = self.prepare(Path(temp))
        self.assertEqual(clinical.pcsm_event.tolist(), [1, 0, 0])
        self.assertTrue(clinical.pcsm_months.isna().all())
        self.assertFalse(clinical.pcsm_eligible.any())
        self.assertEqual(clinical.source_followup_months_from_rrp.tolist(), [60, 72, 84])
        # PSA control with a later systemic interval must not be misclassified as a permanent non-event.
        self.assertEqual(clinical.source_systemic_interval_months_from_psa.iloc[1], 120)
        self.assertNotIn('systemic_event', clinical)
        self.assertNotIn('bcr_event', clinical)

    def test_incomplete_duplicate_or_clinically_conflicting_panels_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cache, _ = self.prepare(root)
            panels = {p: read_matrix(cache.root / f'{nakagawa.SERIES}-{p}_series_matrix.txt.gz',
                                    series=nakagawa.SERIES, platform=p)[:2] for p in nakagawa.PANELS}
            meta, expression = panels['GPL5873']
            for kind in ('missing', 'duplicate', 'conflict'):
                bad = meta.copy()
                values = expression.copy()
                if kind == 'missing':
                    bad = bad.iloc[:-1]
                    values = values.iloc[:, :-1]
                elif kind == 'duplicate':
                    bad.loc[1, 'title'] = bad.title.iloc[0]
                else:
                    bad.loc[1, 'geo__prostate_cancer_specific_death'] = 'No'
                with self.subTest(kind=kind), self.assertRaises(ValueError):
                    nakagawa.join_panels({**panels, 'GPL5873': (bad, values)})

    def test_full_export_and_byte_identical_offline_rebuild(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cache, _ = self.prepare(root)
            roots = [root / p for p in ('processed', 'interim', 'anndata')]
            report = nakagawa_pipeline.build(cache, *roots)
            self.assertEqual(report['assays'], 6)
            self.assertEqual(report['counts']['primary_profiles'], 3)
            self.assertEqual(report['counts']['source_coded_pcsm_events'], 1)
            self.assertEqual(check_dataset(roots[0] / nakagawa.SNAPSHOT, roots[2] / nakagawa.SNAPSHOT)['errors'], [])
            copies = [root / ('offline-' + p) for p in ('processed', 'interim', 'anndata')]
            with patch('urllib.request.urlopen', side_effect=AssertionError('offline HTTP')):
                nakagawa_pipeline.build(SourceCache(cache.root, offline=True, accessions=nakagawa.ACCESSIONS), *copies)
            for a, b in zip(roots, copies):
                for path in (a / nakagawa.SNAPSHOT).iterdir():
                    self.assertEqual(sha256(path), sha256(b / nakagawa.SNAPSHOT / path.name))


if __name__ == '__main__':
    unittest.main()
