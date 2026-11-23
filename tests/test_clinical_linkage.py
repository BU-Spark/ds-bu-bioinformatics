import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from ds_bu_bioinformatics import sun
from ds_bu_bioinformatics.clinical_linkage import apply
from ds_bu_bioinformatics.geo_cohort_pipeline import build
from ds_bu_bioinformatics.output_check import check_dataset
from ds_bu_bioinformatics.sources import SourceCache, sha256
from test_belfast_sun import inputs, transport


def supplement(root):
    path = root / 'clinical.csv'
    path.write_text('sample_id,endpoint,event,time,time_unit\nGSM1,recurrence,1,1.5,years\nGSM2,recurrence,0,24,months\n')
    provenance = root / 'provenance.json'
    provenance.write_text(json.dumps({'study_accession': sun.SERIES, 'source': 'Synthetic clinical linkage fixture',
        'endpoint_definitions': {'recurrence': {'definition': 'Deposited recurrent/non-recurrent class; not relabeled BCR or DFS', 'time_origin': 'surgery',
            'censoring': 'last disease-free follow-up', 'sample_linkage': 'exact GSM identifiers in fixture'}}}))
    return path, provenance


class ClinicalLinkageTest(unittest.TestCase):
    def test_matching_numbers_do_not_allow_different_time_origin(self):
        clinical, endpoints, *_ = apply(self.clinical, sun.ENDPOINTS, self.records, self.provenance,
                                        accession=sun.SERIES, source_definitions=self.definitions)
        definitions = {p: {**d, 'time_origin': 'surgery', 'censoring': 'last disease-free follow-up'}
                       for p, d in self.definitions.items()}
        provenance = json.loads(self.provenance.read_text())
        provenance['endpoint_definitions']['recurrence']['time_origin'] = 'start of radiotherapy'
        self.provenance.write_text(json.dumps(provenance))
        with self.assertRaisesRegex(ValueError, 'incompatible clinical time_origin'):
            apply(clinical, endpoints, self.records, self.provenance,
                  accession=sun.SERIES, source_definitions=definitions)

    def test_rejected_supplement_does_not_poison_corrected_retry(self):
        original = self.provenance.read_text()
        provenance = json.loads(original)
        provenance['endpoint_definitions']['recurrence']['definition'] = 'a different endpoint'
        self.provenance.write_text(json.dumps(provenance))
        roots = [self.root / p for p in ('processed', 'interim', 'anndata')]
        with self.assertRaisesRegex(ValueError, 'incompatible clinical definition'):
            build(self.cache, *roots, cohort=sun, clinical_file=self.records, clinical_provenance=self.provenance)
        self.assertFalse(any(name.startswith('clinical-linkage/') for name in self.cache.manifest['sources']))
        self.provenance.write_text(original)
        report = build(self.cache, *roots, cohort=sun, clinical_file=self.records, clinical_provenance=self.provenance)
        self.assertEqual(report['complete_clinical_records']['recurrence'], 2)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cache = SourceCache(self.root / 'raw' / sun.SNAPSHOT, accessions=[sun.SERIES, sun.PLATFORM])
        with transport(inputs(sun)), contextlib.redirect_stdout(io.StringIO()):
            prepared = sun.prepare(self.cache)
            self.clinical = prepared[0]
            self.definitions = prepared[-1]["endpoints"]
        self.records, self.provenance = supplement(self.root)

    def test_partial_exact_linkage_converts_years_preserves_unknown_subject_and_revalidates(self):
        clinical, endpoints, audit, _ = apply(self.clinical, sun.ENDPOINTS, self.records, self.provenance, accession=sun.SERIES, source_definitions=self.definitions)
        self.assertEqual(clinical.recurrence_months.iloc[:2].tolist(), [18, 24])
        self.assertEqual(clinical.recurrence_eligible.tolist(), [True, True, False])
        self.assertTrue(pd.isna(clinical.recurrence_event.iloc[2]))
        self.assertEqual(audit.time_unit.tolist(), ['years', 'months'])
        self.assertEqual(endpoints['recurrence'](clinical.iloc[0].to_dict()), (1, 18, 'eligible'))
        self.assertEqual(self.clinical.geo__recurrence_status.tolist(), clinical.geo__recurrence_status.tolist())

    def test_conflicting_status_unknown_subject_duplicate_or_undeclared_endpoint_fail(self):
        original = self.records.read_text()
        for changed in [original.replace('GSM1,recurrence,1', 'GSM1,recurrence,0'),
                        original.replace('GSM1,', 'GSM404,'),
                        original + original.splitlines()[1] + '\n',
                        original.replace(',recurrence,', ',os,'),
                        original.replace(',1.5,years', ',1.5,days'),
                        original.replace(',1.5,years', ',,years'),
                        original.replace(',1.5,years', ',-1,years'),
                        original.replace(',1.5,years', ',inf,years')]:
            self.records.write_text(changed)
            with self.subTest(payload=changed), self.assertRaises(ValueError):
                apply(self.clinical, sun.ENDPOINTS, self.records, self.provenance, accession=sun.SERIES, source_definitions=self.definitions)

    def test_known_duration_cannot_be_silently_replaced(self):
        clinical, endpoints, *_ = apply(self.clinical, sun.ENDPOINTS, self.records, self.provenance, accession=sun.SERIES, source_definitions=self.definitions)
        self.definitions['recurrence'].update(time_origin='surgery', censoring='last disease-free follow-up')
        self.records.write_text(self.records.read_text().replace('1.5,years', '2,years'))
        with self.assertRaisesRegex(ValueError, 'duration conflicts'):
            apply(clinical, endpoints, self.records, self.provenance, accession=sun.SERIES, source_definitions=self.definitions)

    def test_wrong_study_or_undocumented_time_origin_fails(self):
        original = json.loads(self.provenance.read_text())
        for changed in [{**original, 'study_accession': 'GSE_OTHER'},
                        {**original, 'endpoint_definitions': {'recurrence': {'definition': 'recurrence'}}}]:
            self.provenance.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):
                apply(self.clinical, sun.ENDPOINTS, self.records, self.provenance, accession=sun.SERIES, source_definitions=self.definitions)

    def test_linked_export_and_cached_offline_replay_survive_original_file_removal(self):
        roots = [self.root / p for p in ('processed', 'interim', 'anndata')]
        report = build(self.cache, *roots, cohort=sun, clinical_file=self.records, clinical_provenance=self.provenance)
        self.assertEqual(report['complete_clinical_records']['recurrence'], 2)
        self.assertFalse(report['analysis_policy']['ordinary_cox_allowed'])
        self.assertEqual(report['counts']['eligible_recurrence_profiles'], 2)
        metadata = json.loads((roots[0] / sun.SNAPSHOT / 'metadata.json').read_text())
        self.assertEqual(metadata['clinical_linkage']['endpoint_definitions']['recurrence']['time_origin'], 'surgery')
        self.assertFalse(metadata['endpoints']['recurrence']['duration_available'])
        self.assertEqual(check_dataset(roots[0] / sun.SNAPSHOT, roots[2] / sun.SNAPSHOT)['errors'], [])
        self.records.unlink()
        self.provenance.unlink()
        copies = [self.root / ('offline-' + p) for p in ('processed', 'interim', 'anndata')]
        with patch('urllib.request.urlopen', side_effect=AssertionError('offline HTTP')):
            build(SourceCache(self.cache.root, offline=True, accessions=[sun.SERIES, sun.PLATFORM]), *copies, cohort=sun)
        for source, copy in zip(roots, copies):
            for path in (source / sun.SNAPSHOT).iterdir():
                self.assertEqual(sha256(path), sha256(copy / sun.SNAPSHOT / path.name))
        clinical = pd.read_csv(roots[0] / sun.SNAPSHOT / 'clinical.csv', keep_default_na=False)
        self.assertEqual(clinical.geo__recurrence_status.tolist(), ['Recurrent', 'Non-Recurrent', 'NA'])

    def test_cache_rejects_changed_local_files_and_checksums(self):
        cached = self.cache.ingest('clinical.csv', self.records, source='fixture')
        self.records.write_text(self.records.read_text().replace('1.5', '2.5'))
        with self.assertRaisesRegex(ValueError, 'Local input changed'):
            self.cache.ingest('clinical.csv', self.records, source='fixture')
        cached.write_text('corrupt')
        with self.assertRaisesRegex(ValueError, 'corrupt cached'):
            self.cache.ingest('clinical.csv', None, source='fixture')
        with self.assertRaises(ValueError):
            self.cache.ingest('../escape', self.records, source='fixture')


if __name__ == '__main__':
    unittest.main()
