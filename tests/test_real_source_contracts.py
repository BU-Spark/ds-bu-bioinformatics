"""Regressions for schemas encountered in the full deposited source files."""
import contextlib
import concurrent.futures
import io
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from ds_bu_bioinformatics import cancermap, cpc_gene, moscow, sun
from ds_bu_bioinformatics.geo_io import read_geo, read_tsv, tsv_feature_inventory
from ds_bu_bioinformatics.exon_arrays import ERHO
from ds_bu_bioinformatics.geo_cohort_pipeline import finite_features
from ds_bu_bioinformatics.sources import SourceCache
from test_belfast_sun import inputs as simple_inputs, transport as simple_transport
from test_remaining_native import gz, matrix, transport, inputs
from scripts.build_with_cached_annotations import ExactAnnotationCache
from ds_bu_bioinformatics.sources import MYGENE_URL


class DepositedSourceContracts(unittest.TestCase):
    def test_annotation_reuse_requires_identical_request_body_and_verified_response(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            peer=SourceCache(root/'peer',accessions=['peer'])
            body=b'q=NM_1&scopes=refseq'
            with transport({}):
                original=peer.fetch('mygene/batch_0000.json',MYGENE_URL,body)
            cache=ExactAnnotationCache(root/'target',offline=True,accessions=['target'])
            with cache.transaction():
                saved=cache.fetch('mygene/batch_0000.json',MYGENE_URL,body)
            self.assertEqual(saved.read_bytes(),original.read_bytes())
            self.assertEqual(cache.manifest['sources']['mygene/batch_0000.json']['retrieved_at_utc'],peer.manifest['sources']['mygene/batch_0000.json']['retrieved_at_utc'])
            changed=ExactAnnotationCache(root/'changed',offline=True,accessions=['changed'])
            with self.assertRaisesRegex(ValueError,'Offline input absent'):
                changed.fetch('mygene/batch_0000.json',MYGENE_URL,b'q=NM_2&scopes=refseq')

    def test_parallel_source_downloads_publish_a_complete_replayable_manifest(self):
        payloads = {f'https://example.org/input-{i}': f'input {i}'.encode() for i in range(20)}
        with tempfile.TemporaryDirectory() as temp, transport(payloads):
            cache = SourceCache(Path(temp), accessions=['test'])
            with cache.transaction(), concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
                paths = list(pool.map(lambda i: cache.fetch(f'input-{i}.txt', f'https://example.org/input-{i}'), range(20)))
            manifest = json.loads(cache.manifest_path.read_text())
            self.assertEqual(set(manifest['sources']), {p.name for p in paths})
            replay = SourceCache(cache.root, offline=True, accessions=['test'])
            replay.verify_all()
            for i in range(20):
                self.assertEqual(replay.fetch(f'input-{i}.txt', f'https://example.org/input-{i}').read_bytes(), payloads[f'https://example.org/input-{i}'])

    def test_erho_supplement_cel_suffix_links_exactly_to_source_description(self):
        import gzip
        payloads = inputs(ERHO)
        payloads[ERHO.SUPPLEMENT_URL] = gz(gzip.decompress(payloads[ERHO.SUPPLEMENT_URL]).decode()
            .replace('GBX.DISC.PCA2\tGBX.DISC.PCA1', 'GBX.DISC.PCA2.CEL\tGBX.DISC.PCA1.CEL'))
        with tempfile.TemporaryDirectory() as temp, transport(payloads):
            prepared = ERHO.prepare(SourceCache(Path(temp), accessions=[ERHO.SERIES, ERHO.PLATFORM]))
        self.assertEqual(prepared[3].columns.tolist(), ['GSM1', 'GSM2'])
        self.assertEqual(prepared[3].iloc[0].tolist(), [2, 4])

    def test_feature_inventory_and_filtered_values_preserve_source_order_and_reject_duplicates(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'matrix.tsv'
            path.write_text('ID_REF\tS2\tS1\na\t1\t2\nb\tNA\t4\nc\t5\t6\n')
            self.assertEqual(tsv_feature_inventory(path).index.tolist(), ['a', 'b', 'c'])
            values = read_tsv(path, feature_filter={'c', 'a'})
            self.assertEqual(values.index.tolist(), ['a', 'c'])
            self.assertEqual(values.columns.tolist(), ['S2', 'S1'])
            path.write_text(path.read_text() + 'b\t7\t8\n')
            with self.assertRaises(ValueError):
                tsv_feature_inventory(path)
            with self.assertRaises(ValueError):
                read_tsv(path, feature_filter={'a'})

    def test_filtered_numeric_loading_keeps_every_annotation_decision_and_unknown_missingness(self):
        source = pd.DataFrame([[1.,2.],[3.,4.]], index=pd.Index(['a','c'],name='feature_id'),columns=['S1','S2'])
        annotations = pd.DataFrame({'feature_id':['a','b','c'], 'gene_id':['1','','2'],
            'gene_symbol':['A','','C'],'mapping_status':['mapped','unmapped','mapped'],
            'annotation_eligible':[True,False,True]})
        expression, genes, mapping = finite_features(source,annotations)
        self.assertEqual(expression.index.tolist(), ['a','c'])
        self.assertEqual(mapping.feature_id.tolist(), ['a','b','c'])
        self.assertTrue(pd.isna(mapping.missing_selected_values.iloc[1]))
        self.assertEqual(mapping.expression_status.iloc[1], 'unmaterialized_unannotated_source')
        with self.assertRaises(ValueError):
            finite_features(source.iloc[:1],annotations)

    def test_moscow_u_namespace_is_not_merged_with_numeric_patient(self):
        meta=pd.DataFrame([{'sample_id':f'GSM{i}','geo_accession':f'GSM{i}',
            'title':title,'platform_id':moscow.PLATFORM,'geo__tissue':'primary tumor',
            'geo__biochemical_recurrence':'brf'} for i,title in enumerate(['1T RNA','U1T RNA'],1)])
        _,samples,_=moscow.clinical_tables(meta,{'GSM1','GSM2'})
        self.assertEqual(samples.patient_id.tolist(),['GSE229904:1','GSE229904:U1'])

    def test_sun_has_both_pg_and_pr_subject_labels(self):
        payloads = simple_inputs(sun)
        import gzip
        payloads[sun.MATRIX_URL] = gz(gzip.decompress(payloads[sun.MATRIX_URL]).decode().replace('tumor PR1', 'tumor PG8'))
        with tempfile.TemporaryDirectory() as temp, simple_transport(payloads):
            prepared = sun.prepare(SourceCache(Path(temp), accessions=[sun.SERIES, sun.PLATFORM]))
        self.assertEqual(prepared[1].patient_id.tolist(), ['PG8', 'PR2', 'PR3'])

    def test_cancermap_low_tumor_and_stroma_are_audited_and_tumor_is_preferred(self):
        meta = pd.DataFrame([{'sample_id': f'GSM{i}', 'geo_accession': f'GSM{i}',
            'platform_id': cancermap.PLATFORM, 'title': f'region {i}', 'organism': 'Homo sapiens',
            'source_name': 'Prostatectomy - ' + material, 'geo__tissue': 'Prostate',
            'geo__material_type': material, 'geo__donor_id': donor}
            for i, material, donor in [(1,'Low Tumor','CI_39'),(2,'Tumor','CI_39'),
                                      (3,'Low Tumor','CI_40'),(4,'Stroma','CI_39')]])
        _, samples, audit = cancermap.clinical_tables(meta)
        self.assertEqual(samples.sample_id.tolist(), ['GSM2', 'GSM3'])
        self.assertEqual(audit.selection_reason.iloc[-1], 'stroma')

    def test_cpc_gene_annotation_columns_are_not_expression_assays(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'processed.tsv'
            path.write_text('GeneID\tSymbol_UCSC\tName_UCSC\tChr_UCSC\tStart_UCSC\tEnd_UCSC\tRefSeq_UCSC\tCPCG0099\tCPCG0100\n'
                            '1\tA1BG\tdescription\tchr19\t58859831\t58874214\tNM_130786\t6.6\t6.3\n')
            expression, annotation = cpc_gene.processed_matrix(path)
            self.assertEqual(expression.columns.tolist(), ['CPCG0099','CPCG0100'])
            self.assertEqual(expression.loc['1'].tolist(), [6.6,6.3])
            self.assertEqual(annotation.Symbol_UCSC.tolist(), ['A1BG'])

    def test_moscow_unregistered_normal_columns_and_rna2_repeats_do_not_become_patients(self):
        records = [{'sample_id':'GSM1','title':'102T RNA','tissue':'primary tumor','biochemical recurrence':'bcr'},
                   {'sample_id':'GSM2','title':'102T RNA2','tissue':'primary tumor','biochemical recurrence':'bcr'},
                   {'sample_id':'GSM3','title':'102N miRNA','tissue':'normal prostate tissue','biochemical recurrence':'NA'},
                   {'sample_id':'GSM4','title':'103T RNA','tissue':'primary tumor','biochemical recurrence':'brf'}]
        payloads = {moscow.MATRIX_URL:matrix(moscow.SERIES,moscow.PLATFORM,records,[]),
                    moscow.COUNTS_URL:gz('gene\t102N\t102T\t102T2\t103T\nENSG000001\t1\t2\t3\t4\n')}
        with tempfile.TemporaryDirectory() as temp, transport(payloads):
            prepared = moscow.prepare(SourceCache(Path(temp),accessions=[moscow.SERIES,moscow.PLATFORM]))
        self.assertEqual(prepared[1].sample_id.tolist(), ['GSM1','GSM4'])
        self.assertEqual(prepared[2].selection_reason.iloc[1], 'repeat_mRNA_same_patient')
        linkage=prepared[5]['expression_linkage.csv']
        excluded=linkage.loc[linkage.source_expression_id.eq('102N')].iloc[0]
        self.assertEqual(excluded.sample_id, '')
        self.assertEqual(excluded.linkage_status, 'excluded_normal_without_GEO_mRNA_record')

    def test_multiline_geo_description_is_preserved_without_weakening_identity_checks(self):
        records=[{'sample_id':'GSM1','title':'source1'}]
        payload=matrix('GSE1','GPL1',records,[['p1',2]],
                       **{'!Sample_description':['first']})
        import gzip
        text=gzip.decompress(payload).decode().replace('!series_matrix_table_begin',
            '!Sample_description\tsecond\n!series_matrix_table_begin')
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'source.gz';path.write_bytes(gz(text))
            meta,_,_=read_geo(path,series='GSE1',platforms={'GPL1'})
        self.assertEqual(meta.source_description.tolist(),['first'])
        self.assertEqual(meta.source_description_2.tolist(),['second'])


if __name__ == '__main__':
    unittest.main()
