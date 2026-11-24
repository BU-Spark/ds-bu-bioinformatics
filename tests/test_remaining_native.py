"""Source-shaped native adapter tests; not a substitute for complete real-source runs."""
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

from ds_bu_bioinformatics import cancermap, cpc_gene, moscow
from ds_bu_bioinformatics import cancermap_pipeline, cpc_gene_pipeline, moscow_pipeline
from ds_bu_bioinformatics import erho_pipeline, karnes_pipeline, den_pipeline, klein_pipeline
from ds_bu_bioinformatics.exon_arrays import ERHO, KARNES, DEN, KLEIN
from ds_bu_bioinformatics.geo_cohort_pipeline import build
from ds_bu_bioinformatics.geo_io import read_geo, read_tsv
from ds_bu_bioinformatics.output_check import check_dataset
from ds_bu_bioinformatics.sources import MYGENE_URL, SourceCache, sha256
from test_mortensen_long import matrix_bytes

NATIVE = [(ERHO, erho_pipeline), (KARNES, karnes_pipeline), (DEN, den_pipeline), (KLEIN, klein_pipeline),
          (cancermap, cancermap_pipeline), (cpc_gene, cpc_gene_pipeline), (moscow, moscow_pipeline)]


def gz(text):
    return gzip.compress(text.encode(), mtime=0)


def matrix(series, platform, records, rows, **extra_headers):
    text = matrix_bytes(series, platform, records, rows).decode()
    extra = ''.join(key + '\t' + '\t'.join(values) + '\n' for key, values in extra_headers.items())
    return gz(text.replace('!series_matrix_table_begin', extra + '!series_matrix_table_begin'))


def platform_payload(platform):
    return ('\n'.join([f'^PLATFORM = {platform}', '!platform_table_begin', 'ID\tgene_assignment',
        'p1\tNM_1 // A /// NR_1 // A', 'p2\tNM_1 // A', 'p3\tNM_2 // B',
        'p4\tNM_1 // A /// NM_2 // B', '!platform_table_end']) + '\n').encode()


def inputs(cohort):
    p = {}
    if cohort in (ERHO, KARNES, DEN, KLEIN):
        records = [{'sample_id': 'GSM1', 'title': 'prostate_cancer_primary_1', **({'metastatic event': '1'} if cohort is ERHO else {})},
                   {'sample_id': 'GSM2', 'title': 'prostate_cancer_primary_2', **({'metastatic event': '0'} if cohort is ERHO else {})}]
        rows = [['p1', 2, 4], ['p2', 5, 8], ['p3', 7, 'NA'], ['p4', 1, 3]]
        p[cohort.MATRIX_URL] = matrix(cohort.SERIES, cohort.PLATFORM, records, [] if cohort is ERHO else rows,
            **{'!Sample_description': ['GBX.DISC.PCA1', 'GBX.DISC.PCA2']})
        if cohort is ERHO:
            p[cohort.SUPPLEMENT_URL] = gz('ID_REF\tGBX.DISC.PCA2\tGBX.DISC.PCA1\np1\t4\t2\np2\t8\t5\np3\tNA\t7\np4\t3\t1\n')
        url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={cohort.PLATFORM}&targ=self&form=text&view=full'
        p[url] = platform_payload(cohort.PLATFORM)
    elif cohort is cancermap:
        records = [{'sample_id': f'GSM{i}', 'title': f'CancerMap region {i}', 'tissue': 'Prostate', 'material_type': tissue, 'donor': donor}
            for i, tissue, donor in [(1, 'malignant', '001'), (2, 'malignant', '001'), (3, 'malignant', '002'),
                                    (4, 'benign', '001'), (5, 'fibroblasts', 'mouse')]]
        p[cohort.MATRIX_URL] = matrix(cohort.SERIES, cohort.PLATFORM, records,
            [['p1', 2, 4, 6, 8, 10], ['p2', 4, 7, 9, 11, 13], ['p3', 1, 3, 'NA', 5, 7], ['p4', 1, 2, 3, 4, 5]],
            **{'!Sample_organism_ch1': ['Homo sapiens'] * 4 + ['Mus musculus']})
        p[f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={cohort.PLATFORM}&targ=self&form=text&view=full'] = platform_payload(cohort.PLATFORM)
    elif cohort is cpc_gene:
        for i, platform in enumerate(cohort.PANELS, 1):
            records = [{'sample_id': f'GSM{i}', 'title': f'prostate tumour CPCG000{i}', 'tissue': 'prostate tumour', 'batch': f'Batch{i}'}]
            p[cohort.MATRIX_URLS[platform]] = matrix(cohort.SERIES, platform, records, [])
        p[cohort.EXPRESSION_URL] = gz('Entrez\tCPCG0002\tCPCG0001\n10_at\t4\t2\n10\t8\t5\n20\tNA\t7\n')
    else:
        records = [{'sample_id': 'GSM1', 'title': '1T RNA', 'tissue': 'primary tumor', 'biochemical recurrence': 'bcr'},
            {'sample_id': 'GSM2', 'title': '2T RNA', 'tissue': 'primary tumor', 'biochemical recurrence': 'brf'},
            {'sample_id': 'GSM3', 'title': '1N RNA', 'tissue': 'normal prostate tissue', 'biochemical recurrence': 'NA'},
            {'sample_id': 'GSM4', 'title': '1T miRNA', 'tissue': 'primary tumor', 'biochemical recurrence': 'bcr'}]
        p[cohort.MATRIX_URL] = matrix(cohort.SERIES, cohort.PLATFORM, records, [])
        p[cohort.COUNTS_URL] = gz('gene_id\t2T\t1N\t1T\nENSG000001.1\t4\t5\t2\nENSG000001.2\t8\t11\t5\nENSG000002\t1\t3\t7\n')
    return p


def transport(payloads):
    class Response(io.BytesIO):
        headers = {}
        def __init__(self, value, url):
            super().__init__(value)
            self.url = url
    def fetch(request, **kwargs):
        if request.full_url == MYGENE_URL:
            from urllib.parse import parse_qs
            args = parse_qs(request.data.decode())
            hits = []
            for q in args['q'][0].split(','):
                number = 1 if q in {'NM_1', 'NR_1', '10', 'ENSG000001'} else 2
                hits.append({'query': q, 'entrezgene': number * 10, 'symbol': 'A' if number == 1 else 'B',
                    'taxid': 9606, 'refseq': {'rna': ['NM_' + str(number), 'NR_' + str(number)]}, 'ensembl': {'gene': f'ENSG00000{number}'}})
            value = json.dumps(hits).encode()
        else:
            value = payloads[request.full_url]
        return Response(value, request.full_url)
    return patch('urllib.request.urlopen', side_effect=fetch)


class NativeCohortsTest(unittest.TestCase):
    def test_documented_cancermap_organ_and_material_are_separate(self):
        record = {'sample_id': 'GSM2482467', 'geo_accession': 'GSM2482467',
                  'platform_id': 'GPL5175', 'title': 'CancerMap tumor region',
                  'organism': 'Homo sapiens', 'source_name': 'primary prostate tumor',
                  'geo__tissue': 'Prostate', 'geo__material_type': 'Tumor', 'geo__donor_id': 'CI_39'}
        xenograft = {**record, 'sample_id': 'GSM9999999', 'geo_accession': 'GSM9999999',
                     'source_name': 'prostate cancer xenograft'}
        _, samples, audit = cancermap.clinical_tables(pd.DataFrame([record, xenograft]))
        self.assertEqual(samples.sample_id.tolist(), ['GSM2482467'])
        self.assertEqual(samples.patient_id.tolist(), ['GSE94767:CI_39'])
        self.assertEqual(audit.selection_reason.iloc[1], 'nonhuman_or_xenograft')

    def prepare(self, root, cohort):
        accessions = getattr(cohort, 'ACCESSIONS', [cohort.SERIES, cohort.PLATFORM])
        cache = SourceCache(root / 'raw' / cohort.SNAPSHOT, accessions=accessions)
        with transport(inputs(cohort)), contextlib.redirect_stdout(io.StringIO()):
            prepared = cohort.prepare(cache)
        return cache, prepared

    def test_every_public_native_pipeline_exports_and_rebuilds_offline_byte_identically(self):
        for cohort, entry in NATIVE:
            with self.subTest(cohort=cohort.SERIES), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                cache, prepared = self.prepare(root, cohort)
                roots = [root / name for name in ('processed', 'interim', 'anndata')]
                report = entry.build(cache, *roots)
                self.assertEqual(report['counts']['primary_profiles'], 2)
                self.assertTrue(report['anndata_round_trip_valid'])
                self.assertEqual(check_dataset(roots[0] / cohort.SNAPSHOT, roots[2] / cohort.SNAPSHOT)['errors'], [])
                saved = pd.read_csv(roots[0] / cohort.SNAPSHOT / 'expression.csv.gz', index_col=0)
                self.assertEqual(saved.iloc[0].tolist(), [2, 4] if cohort is not cancermap else [2, 6])
                self.assertFalse(prepared[0][next(iter(cohort.ENDPOINTS)) + '_eligible'].any())
                copies = [root / ('offline-' + name) for name in ('processed', 'interim', 'anndata')]
                with patch('urllib.request.urlopen', side_effect=AssertionError('offline HTTP')):
                    entry.build(SourceCache(cache.root, offline=True,
                        accessions=getattr(cohort, 'ACCESSIONS', [cohort.SERIES, cohort.PLATFORM])), *copies)
                for a, b in zip(roots, copies):
                    for file in (a / cohort.SNAPSHOT).iterdir():
                        self.assertEqual(sha256(file), sha256(b / cohort.SNAPSHOT / file.name), file.name)

    def test_cancermap_uses_source_donor_and_never_averages_regions(self):
        with tempfile.TemporaryDirectory() as temp:
            _, (clinical, samples, audit, source, *_rest) = self.prepare(Path(temp), cancermap)
        self.assertEqual(samples.sample_id.tolist(), ['GSM1', 'GSM3'])
        self.assertEqual(samples.patient_id.tolist(), ['GSE94767:001', 'GSE94767:002'])
        self.assertEqual(audit.selection_reason.tolist(), ['primary_region_lowest_gsm', 'other_primary_region_same_donor',
            'primary_region_lowest_gsm', 'benign', 'nonhuman_or_xenograft'])
        np.testing.assert_array_equal(source.iloc[0], [2, 6])
        with self.assertRaises(ValueError):
            cancermap.clinical_tables(audit.drop(columns='geo__donor'))

    def test_cpc_exact_processed_matrix_patient_link_and_preserved_duplicate_gene_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            _, (clinical, samples, audit, source, genes, *_rest) = self.prepare(Path(temp), cpc_gene)
        self.assertEqual(samples.patient_id.tolist(), ['CPCG0001', 'CPCG0002'])
        self.assertEqual(source.index.tolist(), ['10_at', '10', '20'])
        self.assertEqual(genes.gene_id.tolist(), ['10', '10', '20'])
        self.assertEqual(clinical.geo__batch.tolist(), ['Batch1', 'Batch2'])

    def test_moscow_excludes_mirna_and_normals_preserves_bcr_and_count_library_sizes(self):
        with tempfile.TemporaryDirectory() as temp:
            _, (clinical, samples, audit, source, genes, extra, processing, details) = self.prepare(Path(temp), moscow)
        self.assertEqual(clinical.bcr_event.tolist(), [1, 0])
        self.assertTrue(clinical.bcr_months.isna().all())
        self.assertEqual(audit.selection_reason.tolist(), ['primary_mRNA', 'primary_mRNA', 'normal_tissue', 'miRNA_assay'])
        self.assertEqual(details['library_sizes'], {'GSM1': 14, 'GSM2': 13})
        self.assertEqual(genes.gene_id.tolist(), ['10', '10', '20'])
        with self.assertRaises(ValueError):
            moscow.outcome({'geo__biochemical_recurrence': 'unknown-code'})

    def test_source_table_duplicate_ids_infinite_values_and_truncation_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for text in ['ID\tS1\tS1\np1\t1\t2\n', 'ID\tS1\np1\t1\np1\t2\n', 'ID\tS1\np1\tinf\n']:
                p = root / 'bad.tsv'
                p.write_text(text)
                with self.assertRaises(ValueError):
                    read_tsv(p)
            payload = inputs(KLEIN)[KLEIN.MATRIX_URL]
            p = root / 'bad.gz'
            p.write_bytes(gzip.compress(gzip.decompress(payload).replace(b'!series_matrix_table_end', b'')))
            with self.assertRaises(ValueError):
                read_geo(p, series=KLEIN.SERIES, platforms={KLEIN.PLATFORM})

    def test_wrong_specimen_and_duplicate_subject_fail_for_each_decipher_cohort(self):
        for cohort in (ERHO, KARNES, DEN, KLEIN):
            with self.subTest(cohort=cohort.SERIES), tempfile.TemporaryDirectory() as temp:
                _, (_, _, audit, *_rest) = self.prepare(Path(temp), cohort)
                for field, value in [('title', 'metastasis_biopsy_1'), ('title', audit.title.iloc[1])]:
                    bad = audit.copy()
                    bad.loc[0, field] = value
                    with self.assertRaises(ValueError):
                        cohort.clinical_tables(bad)


if __name__ == '__main__':
    unittest.main()
