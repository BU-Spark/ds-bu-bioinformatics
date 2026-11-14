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

from ds_bu_bioinformatics import mortensen, long_cohort
from ds_bu_bioinformatics import mortensen_pipeline, long_pipeline
from ds_bu_bioinformatics.cambridge import read_matrix
from ds_bu_bioinformatics.geo_cohort_pipeline import add_outcomes, finite_features, validate
from ds_bu_bioinformatics.sources import SourceCache, sha256


def matrix_bytes(series, platform, records, rows):
    result = [f'!Series_geo_accession\t"{series}"']
    for key, values in [('!Sample_geo_accession', [r['sample_id'] for r in records]),
                        ('!Sample_title', [r['title'] for r in records]),
                        ('!Sample_platform_id', [platform] * len(records)),
                        ('!Sample_data_processing', ['deposited signal'] * len(records))]:
        result.append('\t'.join([key, *values]))
    keys = list(dict.fromkeys(k for r in records for k in r if k not in {'sample_id', 'title'}))
    for key in keys:
        result.append('\t'.join(['!Sample_characteristics_ch1', *[key + ': ' + r[key] if key in r else '' for r in records]]))
    result += ['!series_matrix_table_begin', '\t'.join(['ID_REF', *[r['sample_id'] for r in records]])]
    result.extend('\t'.join(map(str, row)) for row in rows)
    result.append('!series_matrix_table_end')
    return ('\n'.join(result) + '\n').encode()


def mortensen_inputs():
    records = [dict(sample_id='GSM1', title='Prostate cancer cells, PC100', tissue='prostate tumor', bcr='YES', bcr_free_time='7'),
               dict(sample_id='GSM2', title='Prostate cancer cells, PC200', tissue='prostate tumor', bcr='NO', bcr_free_time='12'),
               dict(sample_id='GSM3', title='Prostate cancer cells, PC300', tissue='prostate tumor', bcr='NO', bcr_free_time='NA'),
               dict(sample_id='GSM4', title='Benign prostate glands adjacent to cancer, PC400', tissue='benign prostate glands', bcr='NA', bcr_free_time='NA')]
    rows = [['p1_at', 2, 4, 6, 8], ['p2_at', 8, 10, 12, 14], ['p3_at', 1, 3, '', 7],
            ['p4_at', 2, 3, 4, 5], ['p5_at', 5, 7, 9, ''], ['p6_at', 1, 2, 3, 4],
            ['p7_at', 1, 2, 3, 4], ['p8_at', 2, 3, 4, 5]]
    platform = ['^PLATFORM = GPL570', '!platform_table_begin',
                'ID\tENTREZ_GENE_ID\tGene Symbol\tRefSeq Transcript ID\tSpecies Scientific Name',
                'p1_at\t10\tA\tNM_1\tHomo sapiens', 'p2_at\t10\tA\tNM_2\tHomo sapiens',
                'p3_at\t20\tB\tNM_3\tHomo sapiens', 'p4_at\t10 /// 20\tA /// B\tNM_1 /// NM_3\tHomo sapiens',
                'p5_at\t30\tC\tNM_5\tHomo sapiens', 'p6_at\t40\tD\tNM_6\tMus musculus',
                'p7_at\t50\tE\tNM_7\tHomo sapiens', 'p8_at\t50\tF\tNM_8\tHomo sapiens', '!platform_table_end']
    return {mortensen.MATRIX_URL: gzip.compress(matrix_bytes('GSE46602', 'GPL570', records, rows), mtime=0),
            mortensen.PLATFORM_URL: ('\n'.join(platform) + '\n').encode()}


def long_inputs():
    records = []
    for gsm, cgsr, patient, reads, status, event, followup, title in [
        ('GSM1','CG1','UTPC1','100','1','0','10','Prostate_UTPC1'),
        ('GSM2','CG1-R','UTPC1','200','1','0','10','Prostate_UTPC1-2'),
        ('GSM3','CG2','MCC2','150','0','NA','24','Prostate_MCC2'),
        ('GSM4','CG3','MCC3','250','1','6','10','Prostate_MCC3')]:
        records.append({'sample_id':gsm,'title':title,'cgsr id':cgsr,'external id':patient,'no of uniquely mapped reads':reads,
                        'tissue':'prostate','bcr':status,'months to bcr':event,'months total f/u':followup})
    rows = [['','','CG1','ALIAS1','CG2','CG3'],['External ID','','UTPC1','UTPC1','MCC2','MCC3'],
            ['Months to BCR','','0','0','NA','99'], ['Months total F/U','','10','10','24','88'],
            ['BCR','','1','1','0','0'],['','','','','',''],['Entrez','SYMBOL','','','',''],
            ['10','A','2','4','6','8'],['10','B','8','10','12','14'],
            ['30','C','1','#NUM!','5','7'],['40','D','NA','5','7','9']]
    return {long_cohort.MATRIX_URL: gzip.compress(matrix_bytes('GSE54460','GPL11154',records,[]),mtime=0),
            long_cohort.FPKM_URL: gzip.compress(('\n'.join('\t'.join(row) for row in rows)+'\n').encode(),mtime=0)}


def cached(root, cohort, inputs):
    class Response(io.BytesIO):
        headers = {}
        def __init__(self, value, url):
            super().__init__(value)
            self.url = url
    def transport(request, **kwargs):
        if request.full_url not in inputs:
            raise AssertionError(request.full_url)
        return Response(inputs[request.full_url], request.full_url)
    cache = SourceCache(root/'raw'/'fixture',accessions=[cohort.SERIES,cohort.PLATFORM])
    return cache, patch('urllib.request.urlopen',side_effect=transport)


class MortensenTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cache, transport = cached(self.root, mortensen, mortensen_inputs())
        with transport:
            self.prepared = mortensen.prepare(self.cache)

    def test_tumor_selection_includes_missing_outcomes_and_excludes_benign(self):
        c,s,a,*_ = self.prepared
        self.assertEqual(s.patient_id.tolist(),['PC100','PC200','PC300'])
        self.assertEqual(a.selected.tolist(),[True,True,True,False])
        self.assertEqual(c.bcr_eligible.tolist(),[True,True,False])
        self.assertEqual(c.bcr_months.iloc[:2].tolist(),[7,12])

    def test_annotation_ambiguity_conflicts_missing_rows_and_separate_probes(self):
        c,s,a,x,g,*_=self.prepared
        exp,genes,mapping=finite_features(x,g)
        self.assertEqual(exp.index.tolist(),['p1_at','p2_at','p5_at'])
        self.assertEqual(genes.gene_id.tolist(),['10','10','30'])
        np.testing.assert_array_equal(exp.to_numpy(),[[2,4,6],[8,10,12],[5,7,9]])
        self.assertEqual(mapping.loc[6,'mapping_status'],'gene_symbol_identity_conflict')
        self.assertFalse(mapping.loc[2,'selected'])  # Missing only in ineligible primary.
        self.assertTrue(mapping.loc[4,'selected'])  # Missing only in excluded benign.
        self.assertEqual(validate(exp,c,s,genes,mapping,a,mortensen.outcome)['counts']['bcr_events'],1)

    def test_unknown_tissue_duplicate_subject_and_status_fail(self):
        c,s,a,*_=self.prepared
        for field,value in [('geo__tissue','other'),('title','Prostate cancer cells, PC200')]:
            bad=a.drop(columns=['patient_id','selected','selection_reason']).copy()
            bad.loc[0,field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):
                mortensen.clinical_tables(bad)
        with self.assertRaises(ValueError):mortensen.outcome({'geo__bcr':'MAYBE','geo__bcr_free_time':'4'})
        with self.assertRaises(ValueError):mortensen.outcome({'geo__bcr':'YES','geo__bcr_free_time':'-1'})
        with self.assertRaises(ValueError):mortensen.outcome({'geo__bcr':'YES','geo__bcr_free_time':'nan'})

    def test_zero_and_missing_time_flagged_without_changing_status(self):
        self.assertEqual(mortensen.outcome({'geo__bcr':'YES','geo__bcr_free_time':'0'}),(1,0.,'zero_time_requires_review'))
        self.assertEqual(mortensen.outcome({'geo__bcr':'NO','geo__bcr_free_time':'NA'}),(0,None,'missing_time'))


class LongTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.cache,transport=cached(self.root,long_cohort,long_inputs())
        with transport:self.prepared=long_cohort.prepare(self.cache)

    def test_assay_alias_linkage_and_outcome_independent_quality_selection(self):
        c,s,a,*_=self.prepared
        self.assertEqual(s.sample_id.tolist(),['GSM2','GSM3','GSM4'])
        self.assertEqual(s.patient_id.tolist(),['UTPC1','MCC2','MCC3'])
        self.assertEqual(s.expression_link_method.tolist(),['unique_remaining_external_id','exact_cgsr_id','exact_cgsr_id'])
        self.assertEqual(a.selected.tolist(),[False,True,True,True])
        # Do not substitute the FPKM metadata (status 0, event time 99).
        self.assertEqual(c.bcr_event.tolist(),[1,0,1])
        self.assertEqual(c.bcr_months.tolist(),[0.,24.,6.])
        self.assertEqual(c.bcr_eligible.tolist(),[False,True,True])
        meta=read_matrix(self.cache.root/'GSE54460_series_matrix.txt.gz',series=long_cohort.SERIES,platform=long_cohort.PLATFORM,allow_empty=True)[0]
        _,_,assays,_=long_cohort.read_fpkm(self.cache.root/'GSE54460_FPKM.txt.gz')
        _,shuffled,_=long_cohort.link_and_select(meta.iloc[::-1],assays.iloc[::-1])
        self.assertEqual(set(shuffled.sample_id),{'GSM2','GSM3','GSM4'})
        meta.loc[0,'geo__no_of_uniquely_mapped_reads']='200'
        _,tied,_=long_cohort.link_and_select(meta,assays)
        self.assertIn('GSM1',set(tied.sample_id))  # Stable GSM tie-break.

    def test_duplicate_gene_ids_preserved_and_numeric_error_excluded_without_fill(self):
        c,s,a,x,g,extra,*_=self.prepared
        exp,genes,mapping=finite_features(x,g)
        self.assertEqual(genes.gene_id.tolist(),['10','10','40'])
        self.assertTrue(genes.feature_id.is_unique)
        np.testing.assert_array_equal(exp.to_numpy(),[[4,6,8],[10,12,14],[5,7,9]])
        self.assertEqual(extra['expression_value_audit.csv'].source_token.tolist(),['#NUM!','NA'])
        self.assertEqual(extra['expression_value_audit.csv'].sample_selected.tolist(),[True,False])
        self.assertEqual(mapping.missing_selected_values.tolist(),[0,0,1,0])
        report=validate(exp,c,s,genes,mapping,a,long_cohort.outcome)
        self.assertEqual(report['counts']['source_coded_bcr_events'],2)
        self.assertEqual(report['counts']['bcr_events'],1)
        self.assertEqual(report['counts']['eligible_bcr_profiles'],2)
        comparison=extra['clinical_source_comparison.csv']
        self.assertFalse(comparison.status_agrees.iloc[3])

    def test_conflicting_or_ambiguous_assay_identity_quality_or_patient_status_fails(self):
        meta=read_matrix(self.cache.root/'GSE54460_series_matrix.txt.gz',series=long_cohort.SERIES,platform=long_cohort.PLATFORM,allow_empty=True)[0]
        _,_,assays,_=long_cohort.read_fpkm(self.cache.root/'GSE54460_FPKM.txt.gz')
        for field,value in [('geo__external_id','WRONG'),('geo__no_of_uniquely_mapped_reads','bad'),('geo__bcr','0')]:
            bad=meta.copy();bad.loc[0,field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):long_cohort.link_and_select(bad,assays)
        aliases=assays.copy();aliases.loc[0,'expression_sample_id']='ALIAS2'
        with self.assertRaisesRegex(ValueError,'ambiguous'):long_cohort.link_and_select(meta,aliases)

    def test_invalid_fpkm_tokens_and_duplicate_assays_fail(self):
        payload=long_inputs()[long_cohort.FPKM_URL]
        text=gzip.decompress(payload).decode()
        for old,new in [('#NUM!','garbage'),('#NUM!','inf'),('#NUM!','-3'),('CG1\tALIAS1','CG1\tCG1')]:
            p=self.root/'bad.txt';p.write_text(text.replace(old,new))
            with self.subTest(new=new),self.assertRaises(ValueError):long_cohort.read_fpkm(p)


class ExportTest(unittest.TestCase):
    def test_both_online_exports_and_byte_identical_offline_rebuilds(self):
        for cohort,entry,inputs in [(mortensen,mortensen_pipeline,mortensen_inputs()),(long_cohort,long_pipeline,long_inputs())]:
            with self.subTest(cohort=cohort.SERIES),tempfile.TemporaryDirectory() as temp:
                root=Path(temp);cache,transport=cached(root,cohort,inputs)
                roots=[root/name for name in ['processed','interim','anndata']]
                with transport:report=entry.build(cache,*roots)
                self.assertTrue(report['export_round_trip_valid']);self.assertTrue(report['anndata_round_trip_valid'])
                saved=ad.read_h5ad(roots[2]/'fixture'/'expression.h5ad')
                self.assertEqual(saved.shape,(3,3))
                self.assertEqual(saved.var.gene_id.tolist(),['10','10','30'] if cohort is mortensen else ['10','10','40'])
                copies=[root/('offline-'+n) for n in ['processed','interim','anndata']]
                with patch('urllib.request.urlopen',side_effect=AssertionError('network in offline build')):
                    entry.build(SourceCache(cache.root,offline=True,accessions=[cohort.SERIES,cohort.PLATFORM]),*copies)
                for a,b in zip(roots,copies):
                    for p in (a/'fixture').iterdir():self.assertEqual(sha256(p),sha256(b/'fixture'/p.name),p.name)
                metadata=json.loads(saved.uns['metadata_json']);self.assertEqual(metadata['geo_series'],cohort.SERIES)
                # A conflict in one output format prevents publication of every format.
                blocked=[root/('blocked-'+n) for n in ['processed','interim','anndata']]
                (blocked[2]/'fixture').mkdir(parents=True)
                (blocked[2]/'fixture'/'keep.txt').write_text('old')
                with self.assertRaisesRegex(ValueError,'Existing outputs differ'):entry.build(cache,*blocked)
                self.assertFalse((blocked[0]/'fixture').exists())
                self.assertFalse((blocked[1]/'fixture').exists())
                self.assertEqual((blocked[2]/'fixture'/'keep.txt').read_text(),'old')

    def test_cli_rejects_overlap_and_snapshot_before_network(self):
        for entry in [mortensen_pipeline,long_pipeline]:
            for args in [['--snapshot','../escape'],['--output-root','data/raw'],['--anndata-root','data/processed']]:
                with self.subTest(args=args),patch('urllib.request.urlopen',side_effect=AssertionError('unexpected HTTP')),self.assertRaises(SystemExit):entry.main(args)

    def test_roundtrip_empty_error_audit_is_supported(self):
        inputs=long_inputs()
        inputs[long_cohort.FPKM_URL]=gzip.compress(gzip.decompress(inputs[long_cohort.FPKM_URL]).replace(b'#NUM!',b'3').replace(b'\tNA\t5\t7\t9',b'\t1\t5\t7\t9'))
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);cache,transport=cached(root,long_cohort,inputs)
            with transport:report=long_pipeline.build(cache,root/'out',root/'interim',root/'anndata')
            self.assertEqual(report['expression_invalid_source_cells'],0)


if __name__=='__main__':unittest.main()
