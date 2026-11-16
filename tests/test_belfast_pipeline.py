import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

from ds_bu_bioinformatics import belfast, belfast_pipeline
from ds_bu_bioinformatics.geo_cohort_pipeline import finite_features
from ds_bu_bioinformatics.module_scores import endpoint
from ds_bu_bioinformatics.output_check import check_dataset
from ds_bu_bioinformatics.sources import SourceCache, sha256
from test_mortensen_long import cached, matrix_bytes


def record(i, bcr, bcr_time, met, met_time, tissue='Primary Prostate Cancer', title=None):
    return {'sample_id': f'GSM{i}', 'title': title or f'S1236-0000{i}', 'tissue': tissue, 'patient age (years)': '70',
            'psa (ng/ml)': '12.5', 't-stage': 'T3', 'gleason grade': '8',
            'bcr event (1=yes, 0=no)': bcr, 'follow-up time (bcr, months)': bcr_time,
            'met event (1=yes, 0=no)': met, 'follow-up time (met, months)': met_time}


RECORDS = [record(1, '1', '20', '1', '30'), record(2, '0', '90', '0', '90'),
           record(3, '0', '60', '0', 'NA'), record(4, '1', '0', '0', '40')]
ROWS = [['s1_at', 2, 4, 6, 8],      # sense exonic, gene A
        ['s2_at', 3, 5, 7, 9],      # sense intronic, also gene A: kept as a separate row
        ['s3_at', 1, 2, 3, 4],      # antisense probe that the annotation labels with gene B
        ['s4_at', 1, 2, 3, 4],      # multi-gene
        ['s5_at', 1, 2, 3, 4],      # no genome match
        ['s6_at', 1, 2, 3, 4],      # sense but unannotated
        ['s7_at', 1, 2, 3, 4],      # E and F share Entrez 50: identity conflict
        ['s8_at', 1, 2, 3, 4],
        ['s9_at', 1, 2, '', 4]]     # mapped, but missing in a selected sample
PLATFORM = ['^PLATFORM = GPL25318', '!Platform_title = Almac', '!platform_table_begin',
            'ID\tOrientation\tGene Symbol\tEntrez Gene\tEnsembl Gene ID\tStrand',
            's1_at\tSense (Fully Exonic)\tA\t10\tENSG1\tForward Strand',
            's2_at\tSense (Includes Intronic)\tA\t10\tENSG1\tForward Strand',
            's3_at\tAntiSense\tB\t20\tENSG2\tReverse Strand',
            's4_at\tSense (Fully Exonic)\tC /// D\t30 /// 40\tENSG3\tForward Strand',
            's5_at\tNo Genome Match\t---\t---\t---\t---',
            's6_at\tSense (Fully Exonic)\t---\t---\tENSG6\tForward Strand',
            's7_at\tSense (Fully Exonic)\tE\t50\tENSG7\tForward Strand',
            's8_at\tSense (Fully Exonic)\tF\t50\tENSG8\tForward Strand',
            's9_at\tSense (Fully Exonic)\tG\t70\tENSG9\tForward Strand',
            'unused_at\tAntiSense\tH\t80\tENSG10\tReverse Strand',
            '!platform_table_end']


def inputs(records=RECORDS, platform=PLATFORM):
    return {belfast.MATRIX_URL: gzip.compress(matrix_bytes('GSE116918', 'GPL25318', records, ROWS), mtime=0),
            belfast.PLATFORM_URL: ('\n'.join(platform) + '\n').encode()}


class BelfastTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.cache, transport = cached(self.root, belfast, inputs())
        with transport:
            self.prepared = belfast.prepare(self.cache)

    def test_both_endpoints_decoded_with_eligibility(self):
        clinical, samples, audit, *_ = self.prepared
        self.assertEqual(samples.patient_id.tolist(), ['S1236-00001', 'S1236-00002', 'S1236-00003', 'S1236-00004'])
        self.assertEqual(clinical.metastasis_event.tolist(), [1, 0, 0, 0])
        self.assertEqual(clinical.metastasis_eligibility_reason.tolist(), ['eligible', 'eligible', 'missing_time', 'eligible'])
        self.assertEqual(clinical.bcr_eligibility_reason.tolist(), ['eligible', 'eligible', 'eligible', 'zero_time_requires_review'])
        self.assertEqual(clinical.metastasis_months.iloc[0], 30)
        self.assertTrue(clinical.metastasis_time_unit.eq('months').all())
        self.assertTrue(audit.selected.all())

    def test_only_unambiguous_sense_probes_are_mapped(self):
        clinical, samples, audit, source, annotations, *_ = self.prepared
        status = annotations.set_index('feature_id').mapping_status
        self.assertEqual(status.to_dict(), {'s1_at': 'mapped', 's2_at': 'mapped', 's3_at': 'antisense_probe',
            's4_at': 'multi_gene', 's5_at': 'no_transcript_match', 's6_at': 'unannotated',
            's7_at': 'gene_symbol_identity_conflict', 's8_at': 'gene_symbol_identity_conflict', 's9_at': 'mapped'})
        expression, genes, mapping = finite_features(source, annotations)
        self.assertEqual(expression.index.tolist(), ['s1_at', 's2_at'])     # s9 excluded for missing expression
        self.assertEqual(genes.gene_symbol.tolist(), ['A', 'A'])
        np.testing.assert_array_equal(expression.to_numpy(), [[2, 4, 6, 8], [3, 5, 7, 9]])

    def test_unexpected_inputs_fail(self):
        clinical, samples, audit, *_ = self.prepared
        base = audit.drop(columns=['patient_id', 'selected', 'selection_reason'])
        for field, value in [('geo__tissue', 'Benign'), ('title', 'S1236-00002')]:
            bad = base.copy()
            bad.loc[0, field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                belfast.clinical_tables(bad)
        with self.assertRaises(ValueError):
            belfast.endpoint_outcome('metastasis')({'geo__met_event_1_yes_0_no': '2', 'geo__follow_up_time_met_months': '5'})
        bad_platform = [line.replace('No Genome Match', 'Sideways') for line in PLATFORM]
        cache, transport = cached(self.root / 'bad', belfast, inputs(platform=bad_platform))
        with transport, self.assertRaisesRegex(ValueError, 'orientation'):
            belfast.prepare(cache)


class BelfastExportTest(unittest.TestCase):
    def test_export_passes_shared_checks_and_rebuilds_offline(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cache, transport = cached(root, belfast, inputs())
            roots = [root / name for name in ('processed', 'interim', 'anndata')]
            with transport:
                report = belfast_pipeline.build(cache, *roots)
            self.assertTrue(report['export_round_trip_valid'] and report['anndata_round_trip_valid'])
            self.assertEqual(report['primary_endpoint'], 'metastasis')
            self.assertEqual(report['endpoints']['metastasis']['events'], 1)
            self.assertEqual(report['endpoints']['metastasis']['eligible'], 3)
            self.assertEqual(report['endpoints']['bcr']['eligible'], 3)
            processed = roots[0] / 'fixture'
            metadata = json.loads((processed / 'metadata.json').read_text())
            self.assertEqual(metadata['primary_endpoint'], 'metastasis')
            clinical = pd.read_csv(processed / 'clinical.csv')
            self.assertEqual(endpoint(clinical, metadata, 'belfast-v1')[0], 'metastasis')   # module-scores choice
            result = check_dataset(processed, roots[2] / 'fixture')
            self.assertEqual(result['errors'], [])
            self.assertEqual(set(result['endpoints']), {'bcr', 'metastasis'})
            saved = ad.read_h5ad(roots[2] / 'fixture' / 'expression.h5ad')
            self.assertEqual(saved.shape, (4, 2))
            copies = [root / ('offline-' + n) for n in ('processed', 'interim', 'anndata')]
            with patch('urllib.request.urlopen', side_effect=AssertionError('network in offline build')):
                belfast_pipeline.build(SourceCache(cache.root, offline=True,
                                                   accessions=[belfast.SERIES, belfast.PLATFORM]), *copies)
            for a, b in zip(roots, copies):
                for p in (a / 'fixture').iterdir():
                    self.assertEqual(sha256(p), sha256(b / 'fixture' / p.name), p.name)

    def test_cli_rejects_overlap_and_snapshot_before_network(self):
        for args in [['--snapshot', '../escape'], ['--output-root', 'data/raw'], ['--anndata-root', 'data/processed']]:
            with self.subTest(args=args), patch('urllib.request.urlopen', side_effect=AssertionError('unexpected HTTP')), \
                    self.assertRaises(SystemExit):
                belfast_pipeline.main(args)


if __name__ == '__main__':
    unittest.main()
