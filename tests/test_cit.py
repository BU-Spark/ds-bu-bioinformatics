import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import pandas as pd

from ds_bu_bioinformatics import cit
from ds_bu_bioinformatics.geo_cohort_pipeline import build
from ds_bu_bioinformatics.sources import SourceCache, write_json, sha256


class CITTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.cache = SourceCache(self.root/'raw'/'cit-fixture', accessions=[cit.SERIES,cit.PLATFORM])
        self.names = [f'OCT_{i:03}_HuGene2.CEL' for i in range(1,5)]
        self.header = ['Source Name','Assay Name','Array Data File','Array Design REF',
                       'Description','Characteristics[patient id]',cit.BCR,cit.TIME,'Unit[time unit]',
                       'Factor Value[bcr event]','Factor Value[months to bcr or last news]',
                       'Unit[time unit]','Term Source REF','Term Accession Number']
        self.rows = [[f'PROST{i:03}',name[:-4],name,cit.PLATFORM,label,patient,status,time,unit,status,time,*tail]
                     for i,(name,label,patient,status,time,unit,tail) in enumerate(zip(self.names,
                         ['Tumoral tissue']*3+['Non tumoral tissue'],
                         ['Patient1','Patient2','Patient3','Patient1'],
                         ['1','0','not available','0'],
                         ['12','0','not applicable','24'],
                         ['month','month','','month'],
                         [['month','EFO','UO:0000035'],['month','EFO','UO:0000035'],[],['month','EFO','UO:0000035']]),1)]
        self.sdrf = self.root/'sdrf.txt'
        self.write_sdrf()

    def write_sdrf(self):
        with self.sdrf.open('w',newline='') as f:
            csv.writer(f,delimiter='\t').writerows([self.header,*self.rows])

    def test_tissue_selection_and_clinical_eligibility(self):
        clinical,samples,audit = cit.read_sdrf(self.sdrf)
        self.assertEqual(samples.patient_id.tolist(),['Patient1','Patient2','Patient3'])
        self.assertEqual(audit.selected.tolist(),[True,True,True,False])
        self.assertEqual(clinical.bcr_eligibility_reason.tolist(),['eligible','zero_time_requires_review','missing_status'])
        self.assertTrue(pd.isna(clinical.bcr_months.iloc[2]))

    def test_rejects_factor_conflicts_duplicate_primary_and_wrong_units(self):
        for row,col,value in [(0,9,'0'),(1,5,'Patient1'),(0,8,'year'),(0,2,'../a.CEL')]:
            original = self.rows[row][col]
            self.rows[row][col] = value
            self.write_sdrf()
            with self.subTest(value=value), self.assertRaises(ValueError):
                cit.read_sdrf(self.sdrf)
            self.rows[row][col] = original

    def test_annotation_does_not_average_or_pick_ambiguous_genes(self):
        path = self.root/'annotations.tsv'
        path.write_text('PROBEID\tENTREZID\tSYMBOL\nA\t1\tONE\nB\t1\tONE\nC\t1\tONE\nC\t2\tTWO\nD\t\t\n')
        frame = cit.read_annotations(path,['B','A','C','D'])
        self.assertEqual(frame.feature_id.tolist(),['B','A','C','D'])
        self.assertEqual(frame.annotation_eligible.tolist(),[True,True,False,False])

    def test_export_and_identical_offline_with_mocked_rma(self):
        self.cache.root.mkdir(parents=True)
        write_json(self.cache.manifest_path, self.cache.manifest)
        record = {'accno':cit.SERIES,'files':[{'path':n,'size':1} for n in self.names]}
        paths = {}
        for n in self.names:
            paths[n] = self.root/n
            paths[n].write_bytes(b'x')
        frame = pd.DataFrame([[1.,2.,3.,4.],[9.,8.,7.,6.]],index=pd.Index(['A','B'],name='feature_id'),columns=self.names)
        annotations = pd.DataFrame({'feature_id':['A','B'],'gene_id':['1','1'],'gene_symbol':['ONE','ONE'],
                                    'annotation_eligible':[True,True],'mapping_status':['mapped','mapped']})
        def fetch(name,url):
            return self.sdrf if name == cit.SDRF else paths[name]
        first = [self.root/n for n in ('processed','interim','anndata')]
        second = [self.root/('rebuild-'+n) for n in ('processed','interim','anndata')]
        with patch.object(self.cache,'fetch',side_effect=fetch), patch.object(self.cache,'json',return_value=record), patch.object(cit,'run_rma',return_value=(frame,annotations,['fixture'])) as rma:
            report = build(self.cache,*first,cohort=cit)
            build(self.cache,*second,cohort=cit)
            self.assertEqual([p.name for p in rma.call_args.args[0]],self.names)
        self.assertEqual(report['counts']['primary_profiles'],3)
        saved = ad.read_h5ad(first[2]/'cit-fixture'/'expression.h5ad')
        self.assertEqual(saved.X.tolist(),[[1.,9.],[2.,8.],[3.,7.]])
        for a,b in zip(first,second):
            for f in (a/'cit-fixture').iterdir():
                self.assertEqual(sha256(f),sha256(b/'cit-fixture'/f.name))


if __name__ == '__main__':
    unittest.main()
