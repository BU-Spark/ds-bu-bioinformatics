import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from ds_bu_bioinformatics.authorized_rnaseq import AuthorizedRNASeq, main
from ds_bu_bioinformatics import cpgea_pipeline, dkfz_eopc_pipeline
from ds_bu_bioinformatics.output_check import check_dataset
from ds_bu_bioinformatics.sources import SourceCache, sha256
from test_remaining_native import transport


def sources(root, key):
    expression = root / 'source_counts.tsv'
    expression.write_text('GeneID\teA1\teB1\teA2\teN1\n10\t2\t4\t6\t8\n20\t3\t5\t7\t9\n30\t1\t2\t3\t4\n')
    samples = root / 'clinical_source.tsv'
    samples.write_text('Assay\tDonor\tExpressionKey\tTissue\tRecurrence\tYears\nS1\t001\teA1\ttumor\tYes\t1\nS3\t002\teB1\ttumor\tNo\t3\nS2\t001\teA2\ttumor\tYes\t1.0\nN1\t001\teN1\tnormal\tNA\tNA\n')
    config = {'study_accession': 'EGAD00001004791' if key == 'dkfz-eopc' else 'PRJCA001124',
        'source': 'Synthetic authorized processed research fixture; no actual controlled data',
        'expression_kind': 'raw_counts', 'expression_scale': 'Original source counts',
        'sample_linkage': 'Exact ExpressionKey to source matrix columns, Donor to patient',
        'specimen_definition': 'Tumor specimens retained; normal tissue excluded',
        'sample_columns': {'sample_id': 'Assay', 'patient_id': 'Donor', 'expression_id': 'ExpressionKey', 'tissue': 'Tissue'},
        'tumor_values': ['tumor'], 'excluded_tissue_values': ['normal'], 'replicate_policy': 'lowest_sample_id',
        'primary_endpoint': 'bcr', 'endpoints': {'bcr': {'definition': 'Synthetic source recurrence',
            'time_origin': 'surgery', 'censoring': 'individual last recurrence-free follow-up',
            'event_column': 'Recurrence', 'event_codes': {'Yes': 1, 'No': 0}, 'time_column': 'Years', 'time_unit': 'years'}}}
    path = root / 'source-config.json'
    path.write_text(json.dumps(config))
    return expression, samples, path, config


class AuthorizedStudyTest(unittest.TestCase):
    def test_source_slashes_and_literal_escape_sequences_round_trip_without_key_collisions(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            expression, samples, config_path, config = sources(root, 'cpgea')
            rows = samples.read_text().replace('Years', 'Time/Years').splitlines()
            samples.write_text(rows[0] + '\tTime%2FYears\n' + ''.join(row + '\toriginal\n' for row in rows[1:]))
            config['endpoints']['bcr']['time_column'] = 'Time/Years'
            config_path.write_text(json.dumps(config))
            cohort = AuthorizedRNASeq('cpgea', expression_file=expression, sample_file=samples, source_config=config_path)
            cache = SourceCache(root / 'raw' / cohort.SNAPSHOT, accessions=cohort.ACCESSIONS)
            roots = [root / p for p in ('processed', 'interim', 'anndata')]
            with transport({}), contextlib.redirect_stdout(io.StringIO()):
                cpgea_pipeline.build(cache, *roots, expression_file=expression, sample_file=samples, source_config=config_path)
            self.assertEqual(check_dataset(roots[0] / cohort.SNAPSHOT, roots[2] / cohort.SNAPSHOT)['errors'], [])
            clinical = pd.read_csv(roots[0] / cohort.SNAPSHOT / 'clinical.csv')
            self.assertEqual(clinical.bcr_months.tolist(), [12, 36])
            self.assertEqual(clinical['source__Time%2FYears'].tolist(), [1, 3])
            self.assertEqual(clinical['source__Time%252FYears'].tolist(), ['original', 'original'])
            metadata = json.loads((roots[0] / cohort.SNAPSHOT / 'metadata.json').read_text())
            self.assertEqual(metadata['policies']['source_column_export_names']['Time/Years'], 'source__Time%2FYears')

    def test_both_dedicated_studies_parse_sources_select_regions_and_replay_offline(self):
        for key, entry in [('dkfz-eopc', dkfz_eopc_pipeline), ('cpgea', cpgea_pipeline)]:
            with self.subTest(cohort=key), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                expression, samples, config, _ = sources(root, key)
                cohort = AuthorizedRNASeq(key, expression_file=expression, sample_file=samples, source_config=config)
                cache = SourceCache(root / 'raw' / cohort.SNAPSHOT, accessions=cohort.ACCESSIONS)
                roots = [root / p for p in ('processed', 'interim', 'anndata')]
                with transport({}), contextlib.redirect_stdout(io.StringIO()):
                    report = entry.build(cache, *roots, expression_file=expression, sample_file=samples, source_config=config)
                self.assertEqual(report['counts']['primary_profiles'], 2)
                self.assertEqual(report['counts']['bcr_events'], 1)
                self.assertEqual(report['other_regions_excluded'], 1)
                self.assertEqual(check_dataset(roots[0] / cohort.SNAPSHOT, roots[2] / cohort.SNAPSHOT)['errors'], [])
                clinical = pd.read_csv(roots[0] / cohort.SNAPSHOT / 'clinical.csv', dtype={'patient_id': str})
                self.assertEqual(clinical.patient_id.tolist(), ['001', '002'])
                self.assertEqual(clinical.bcr_months.tolist(), [12, 36])
                self.assertEqual(clinical.source__Years.tolist(), [1, 3])
                for path in (expression, samples, config):
                    path.unlink()
                copies = [root / ('offline-' + p) for p in ('processed', 'interim', 'anndata')]
                with patch('urllib.request.urlopen', side_effect=AssertionError('offline HTTP')):
                    entry.build(SourceCache(cache.root, offline=True, accessions=cohort.ACCESSIONS), *copies)
                for a, b in zip(roots, copies):
                    for path in (a / cohort.SNAPSHOT).iterdir():
                        self.assertEqual(sha256(path), sha256(b / cohort.SNAPSHOT / path.name))

    def test_missing_followup_is_exported_without_inventing_survival(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            expression, samples, config_path, config = sources(root, 'cpgea')
            config['endpoints']['bcr']['time_column'] = None
            config_path.write_text(json.dumps(config))
            cohort = AuthorizedRNASeq('cpgea', expression_file=expression, sample_file=samples, source_config=config_path)
            cache = SourceCache(root / 'raw' / cohort.SNAPSHOT, accessions=cohort.ACCESSIONS)
            with transport({}), contextlib.redirect_stdout(io.StringIO()):
                clinical = cohort.prepare(cache)[0]
            self.assertEqual(clinical.bcr_event.tolist(), [1, 0])
            self.assertTrue(clinical.bcr_months.isna().all())
            self.assertFalse(clinical.bcr_eligible.any())

    def test_wrong_accession_unknown_tissue_duplicate_areas_and_conflicting_outcomes_fail(self):
        for variant in ('accession', 'tissue', 'duplicate', 'conflict', 'unmatched', 'negative_counts'):
            with self.subTest(variant=variant), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                expression, samples, config_path, config = sources(root, 'dkfz-eopc')
                if variant == 'accession':
                    config['study_accession'] = 'PRAD_DIFFERENT_DKFZ_COHORT'
                elif variant == 'tissue':
                    config['excluded_tissue_values'] = []
                elif variant == 'duplicate':
                    config['replicate_policy'] = 'error'
                elif variant == 'conflict':
                    samples.write_text(samples.read_text().replace('S2\t001\teA2\ttumor\tYes', 'S2\t001\teA2\ttumor\tNo'))
                elif variant == 'unmatched':
                    expression.write_text(expression.read_text().replace('eA1', 'unknown-assay'))
                else:
                    expression.write_text(expression.read_text().replace('10\t2\t4', '10\t-2\t4'))
                config_path.write_text(json.dumps(config))
                cohort = AuthorizedRNASeq('dkfz-eopc', expression_file=expression, sample_file=samples, source_config=config_path)
                cache = SourceCache(root / 'raw' / cohort.SNAPSHOT, accessions=cohort.ACCESSIONS)
                with transport({}), contextlib.redirect_stdout(io.StringIO()), self.assertRaises(ValueError):
                    cohort.prepare(cache)

    def test_first_run_requires_real_source_files_and_rejects_overlap(self):
        for key in ('dkfz-eopc', 'cpgea'):
            with tempfile.TemporaryDirectory() as temp:
                cases = [['--data-dir', str(Path(temp) / 'no-inputs')],
                         ['--snapshot', '../escape'], ['--output-root', 'data/raw']]
                for args in cases:
                    with self.subTest(key=key, args=args), self.assertRaises(SystemExit), \
                         contextlib.redirect_stderr(io.StringIO()):
                        main(key, args)


if __name__ == '__main__':
    unittest.main()
