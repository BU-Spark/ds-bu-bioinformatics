"""Synthetic source-shaped tests, not proof of full real-cohort reproducibility."""
import gzip
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

from ds_bu_bioinformatics import belfast, belfast_pipeline, sun, sun_pipeline
from ds_bu_bioinformatics.geo_cohort_pipeline import finite_features, validate
from ds_bu_bioinformatics.output_check import check_dataset
from ds_bu_bioinformatics.sources import SourceCache, sha256


def inputs(cohort):
    if cohort is belfast:
        titles = ['S1236-00001', 'S1236-00002', 'S1236-00003']
        fields = {
            'tissue': ['Primary Prostate Cancer'] * 3,
            'bcr event (1=yes, 0=no)': ['1', '0', '1'],
            'follow-up time (bcr, months)': ['5', '20', 'NA'],
            'met event (1=yes, 0=no)': ['0', '1', '1'],
            'follow-up time (met, months)': ['30', '12', '0']}
        header = ['ID', 'Ensembl Gene ID', 'Gene Symbol', 'Entrez Gene']
        rows = [['p1', 'ENSG10', 'A', '10'], ['p2', 'ENSG10', 'A', '10'],
                ['p3', 'ENSG20', 'B', '20'], ['p4', 'ENSG30 /// ENSG40', 'C /// D', '30 /// 40'],
                ['p5', 'ENSG50', 'E', '50'], ['p6', 'ENSG50', 'F', '50']]
    else:
        titles = [f'Prostate cancer primary tumor PR{i}' for i in range(1, 4)]
        fields = {'tissue': ['Prostate cancer primary tumor'] * 3,
                  'recurrence status': ['Recurrent', 'Non-Recurrent', 'NA']}
        header = ['ID', 'ENTREZ_GENE_ID', 'Gene Symbol', 'Species Scientific Name', 'RefSeq Transcript ID']
        rows = [['p1', '10', 'A', 'Homo sapiens', 'NM_1'], ['p2', '10', 'A', 'Homo sapiens', 'NM_2'],
                ['p3', '20', 'B', 'Homo sapiens', 'NM_3'],
                ['p4', '30 /// 40', 'C /// D', 'Homo sapiens', 'NM_4'],
                ['p5', '50', 'E', 'Homo sapiens', 'NM_5'], ['p6', '50', 'F', 'Homo sapiens', 'NM_6']]
    records = [f'!Series_geo_accession\t{cohort.SERIES}', '!Sample_geo_accession\tGSM1\tGSM2\tGSM3',
               '!Sample_title\t' + '\t'.join(titles),
               '!Sample_platform_id\t' + '\t'.join([cohort.PLATFORM] * 3),
               '!Sample_data_processing\tdeposited\tdeposited\tdeposited']
    records += ['!Sample_characteristics_ch1\t' + '\t'.join(f'{key}: {v}' for v in values)
                for key, values in fields.items()]
    records += ['!series_matrix_table_begin', 'ID_REF\tGSM1\tGSM2\tGSM3',
                'p1\t2\t4\t6', 'p2\t8\t10\t12', 'p3\t1\t3\tNA',
                'p4\t2\t3\t4', 'p5\t3\t4\t5', 'p6\t4\t5\t6', '!series_matrix_table_end']
    platform = [f'^PLATFORM = {cohort.PLATFORM}', '!platform_table_begin', '\t'.join(header),
                *['\t'.join(row) for row in rows], '!platform_table_end']
    return {cohort.MATRIX_URL: gzip.compress(('\n'.join(records) + '\n').encode(), mtime=0),
            cohort.PLATFORM_URL: ('\n'.join(platform) + '\n').encode()}


def transport(payloads):
    class Response(io.BytesIO):
        headers = {}

        def __init__(self, payload, url):
            super().__init__(payload)
            self.url = url

    def fetch(request, **kwargs):
        return Response(payloads[request.full_url], request.full_url)

    return patch('urllib.request.urlopen', side_effect=fetch)


class NewCohortsTest(unittest.TestCase):
    def prepare(self, root, cohort, payloads=None):
        cache = SourceCache(root / 'raw' / cohort.SNAPSHOT, accessions=[cohort.SERIES, cohort.PLATFORM])
        with transport(payloads or inputs(cohort)):
            prepared = cohort.prepare(cache)
        return cache, prepared

    def test_belfast_separate_endpoint_times_missing_and_zero_are_retained(self):
        with tempfile.TemporaryDirectory() as temp:
            _, (clinical, samples, audit, *_rest) = self.prepare(Path(temp), belfast)
        self.assertEqual(clinical.bcr_event.tolist(), [1, 0, 1])
        self.assertEqual(clinical.bcr_months.iloc[:2].tolist(), [5, 20])
        self.assertTrue(pd.isna(clinical.bcr_months.iloc[2]))
        self.assertEqual(clinical.bcr_eligible.tolist(), [True, True, False])
        self.assertEqual(clinical.metastasis_months.tolist(), [30, 12, 0])
        self.assertEqual(clinical.metastasis_eligible.tolist(), [True, True, False])
        self.assertEqual(clinical.metastasis_eligibility_reason.iloc[2], 'zero_time_requires_review')
        self.assertEqual(len(samples), 3)
        self.assertTrue(audit.selected.all())

    def test_sun_status_only_never_invents_censoring_times(self):
        with tempfile.TemporaryDirectory() as temp:
            _, (clinical, samples, *_rest) = self.prepare(Path(temp), sun)
        self.assertEqual(clinical.recurrence_event.iloc[:2].tolist(), [1, 0])
        self.assertTrue(pd.isna(clinical.recurrence_event.iloc[2]))
        self.assertTrue(clinical.recurrence_months.isna().all())
        self.assertFalse(clinical.recurrence_eligible.any())
        self.assertEqual(clinical.recurrence_eligibility_reason.tolist(), ['missing_time', 'missing_time', 'missing_status'])
        self.assertEqual(samples.patient_id.tolist(), ['PR1', 'PR2', 'PR3'])
        self.assertNotIn('bcr_event', clinical)

    def test_feature_audit_preserves_distinct_probes_and_excludes_missing_even_in_ineligible_subject(self):
        for cohort in (belfast, sun):
            with self.subTest(cohort=cohort.SERIES), tempfile.TemporaryDirectory() as temp:
                _, (_, _, _, source, annotations, *_rest) = self.prepare(Path(temp), cohort)
                expression, genes, mapping = finite_features(source, annotations)
                self.assertEqual(genes.feature_id.tolist(), ['p1', 'p2'])
                self.assertEqual(genes.gene_symbol.tolist(), ['A', 'A'])
                np.testing.assert_array_equal(expression.to_numpy(), [[2, 4, 6], [8, 10, 12]])
                self.assertEqual(mapping.mapping_status.iloc[3], 'invalid_or_ambiguous_annotation')
                self.assertEqual(mapping.mapping_status.iloc[4:].tolist(), ['gene_symbol_identity_conflict'] * 2)
                self.assertEqual(mapping.missing_selected_values.iloc[2], 1)

    def test_wrong_subject_duplicate_patient_tissue_or_status_fail(self):
        for cohort in (belfast, sun):
            with self.subTest(cohort=cohort.SERIES), tempfile.TemporaryDirectory() as temp:
                _, (clinical, *_rest) = self.prepare(Path(temp), cohort)
                for field, value in [('title', clinical.title.iloc[1]), ('geo__tissue', 'benign'),
                                     ('sample_id', clinical.sample_id.iloc[1]), ('platform_id', 'GPL_OTHER')]:
                    bad = clinical.copy()
                    bad.loc[0, field] = value
                    with self.subTest(field=field), self.assertRaises(ValueError):
                        cohort.clinical_tables(bad)
                key = 'geo__bcr_event_1_yes_0_no' if cohort is belfast else 'geo__recurrence_status'
                with self.assertRaises(ValueError):
                    cohort.outcome({key: 'MAYBE'})
                with self.assertRaises(ValueError):
                    cohort.clinical_tables(clinical.drop(columns=[key]))

    def test_belfast_negative_nonfinite_or_malformed_time_fails(self):
        for bad in ('-1', 'nan', 'inf', 'bad'):
            with self.subTest(time=bad), self.assertRaises(ValueError):
                belfast.metastasis_outcome({'geo__met_event_1_yes_0_no': '1',
                                           'geo__follow_up_time_met_months': bad})

    def test_secondary_endpoint_tampering_fails_validation(self):
        with tempfile.TemporaryDirectory() as temp:
            _, (clinical, samples, audit, source, annotations, *_rest) = self.prepare(Path(temp), belfast)
            expression, genes, mapping = finite_features(source, annotations)
            clinical.loc[0, 'metastasis_months'] = 999
            with self.assertRaises(AssertionError):
                validate(expression, clinical, samples, genes, mapping, audit, belfast.outcome,
                         endpoints=belfast.ENDPOINTS)

    def test_exports_pass_shared_checker_and_rebuild_byte_identically_offline(self):
        for cohort, entry in ((belfast, belfast_pipeline), (sun, sun_pipeline)):
            with self.subTest(cohort=cohort.SERIES), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                cache = SourceCache(root / 'raw' / cohort.SNAPSHOT, accessions=[cohort.SERIES, cohort.PLATFORM])
                roots = [root / name for name in ('processed', 'interim', 'anndata')]
                with transport(inputs(cohort)):
                    report = entry.build(cache, *roots)
                self.assertTrue(report['anndata_round_trip_valid'])
                checked = check_dataset(roots[0] / cohort.SNAPSHOT, roots[2] / cohort.SNAPSHOT)
                self.assertEqual(checked['errors'], [])
                data = ad.read_h5ad(roots[2] / cohort.SNAPSHOT / 'expression.h5ad')
                self.assertEqual(data.shape, (3, 2))
                metadata = json.loads(data.uns['metadata_json'])
                self.assertEqual(metadata['primary_endpoint'], next(iter(cohort.ENDPOINTS)))
                self.assertEqual(set(metadata['endpoints']), set(cohort.ENDPOINTS))
                self.assertEqual(set(checked['endpoints']), set(cohort.ENDPOINTS))
                copies = [root / ('offline-' + name) for name in ('processed', 'interim', 'anndata')]
                with patch('urllib.request.urlopen', side_effect=AssertionError('offline HTTP')):
                    offline = SourceCache(cache.root, offline=True, accessions=[cohort.SERIES, cohort.PLATFORM])
                    entry.build(offline, *copies)
                for online_root, offline_root in zip(roots, copies):
                    for file in (online_root / cohort.SNAPSHOT).iterdir():
                        self.assertEqual(sha256(file), sha256(offline_root / cohort.SNAPSHOT / file.name), file.name)
                # Existing output conflicts must not partially publish the other formats.
                blocked = [root / ('blocked-' + name) for name in ('processed', 'interim', 'anndata')]
                (blocked[2] / cohort.SNAPSHOT).mkdir(parents=True)
                (blocked[2] / cohort.SNAPSHOT / 'keep.txt').write_text('original')
                with self.assertRaisesRegex(ValueError, 'Existing outputs differ'):
                    entry.build(offline, *blocked)
                self.assertFalse((blocked[0] / cohort.SNAPSHOT).exists())
                self.assertFalse((blocked[1] / cohort.SNAPSHOT).exists())

    def test_source_schema_changes_and_cache_corruption_fail_before_publication(self):
        for cohort, entry in ((belfast, belfast_pipeline), (sun, sun_pipeline)):
            with self.subTest(cohort=cohort.SERIES), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                payloads = inputs(cohort)
                payloads[cohort.PLATFORM_URL] = payloads[cohort.PLATFORM_URL].replace(b'Gene Symbol', b'Changed Symbol')
                with self.assertRaisesRegex(ValueError, 'schema'):
                    self.prepare(root, cohort, payloads)
            with tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                cache, _ = self.prepare(root, cohort)
                file = cache.root / (cohort.SERIES + '_series_matrix.txt.gz')
                file.write_bytes(b'corrupt')
                with self.assertRaisesRegex(ValueError, 'corrupt cached'):
                    entry.build(cache, root / 'out', root / 'interim', root / 'anndata')
                self.assertFalse((root / 'out').exists())

    def test_cli_invalid_paths_fail_without_network(self):
        for entry in (belfast_pipeline, sun_pipeline):
            for args in (['--snapshot', '../escape'], ['--output-root', 'data/raw'],
                         ['--anndata-root', 'data/processed']):
                with self.subTest(entry=entry.__name__, args=args), self.assertRaises(SystemExit), \
                     patch('urllib.request.urlopen', side_effect=AssertionError('unexpected HTTP')):
                    entry.main(args)


if __name__ == '__main__':
    unittest.main()
