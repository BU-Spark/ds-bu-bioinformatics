import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from ds_bu_bioinformatics import dkfz_public as cohort
from ds_bu_bioinformatics import dkfz_public_pipeline as entry
from ds_bu_bioinformatics.output_check import check_dataset
from ds_bu_bioinformatics.sources import SourceCache, sha256
from test_remaining_native import transport


class PublicDKFZTest(unittest.TestCase):
    def test_exact_source_identity_regional_selection_duplicate_feature_rows_and_offline_replay(self):
        payloads={cohort.TEXT_BASE+'meta_study.txt':b'cancer_study_identifier: prostate_dkfz_2018\npmid: 30537516\n',
            cohort.TEXT_BASE+'meta_mrna_seq_rpkm.txt':b'cancer_study_identifier: prostate_dkfz_2018\ndatatype: CONTINUOUS\n',
            cohort.BASE+'data_clinical_patient.txt':b'PATIENT_ID\tAGE_GROUP\tBCR_STATUS\tTIME_FROM_SURGERY_TO_BCR_LASTFU\nP1\tEO\t1\t23.5\nP2\tEO\t0\t14\nP3\tEO\tNA\tNA\n',
            cohort.BASE+'data_clinical_sample.txt':b'PATIENT_ID\tSAMPLE_ID\tSOURCE\tSAMPLE_CLASS\tLOCALIZED_TUMOR\nP1\tS1\tICGC\tTumor\tTRUE\nP1\tS2\tICGC\tTumor\tTRUE\nP2\tS3\tICGC\tTumor\tTRUE\nP3\tS4\tICGC\tNormal\tTRUE\n',
            cohort.BASE+'data_mrna_seq_rpkm.txt':b'Hugo_Symbol\tEntrez_Gene_Id\tS2\tS3\tS1\nB\t\t2\t4\t6\nB\t\t8\t10\t12\nXLOC_1\t\t1\t3\t5\n'}
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);cache=SourceCache(root/'raw'/cohort.SNAPSHOT,accessions=cohort.ACCESSIONS)
            roots=[root/n for n in ['processed','interim','anndata']]
            with transport(payloads),contextlib.redirect_stdout(io.StringIO()):report=entry.build(cache,*roots)
            self.assertEqual(report['counts']['primary_profiles'],2)
            self.assertEqual(report['counts']['retained_feature_rows'],2)
            self.assertEqual(report['other_tumor_regions_excluded'],1)
            self.assertEqual(check_dataset(roots[0]/cohort.SNAPSHOT,roots[2]/cohort.SNAPSHOT)['errors'],[])
            expression=pd.read_csv(roots[0]/cohort.SNAPSHOT/'expression.csv.gz',index_col=0)
            self.assertEqual(expression.index.tolist(),['B_row1','B_row2'])
            self.assertEqual(expression.values.tolist(),[[6.,4.],[12.,10.]])
            clinical=pd.read_csv(roots[0]/cohort.SNAPSHOT/'clinical.csv')
            self.assertEqual(clinical.bcr_event.tolist(),[1,0]);self.assertTrue(clinical.bcr_months.isna().all())
            copies=[root/('replay-'+n) for n in ['processed','interim','anndata']]
            with patch('urllib.request.urlopen',side_effect=AssertionError('Offline HTTP')):
                entry.build(SourceCache(cache.root,offline=True,accessions=cohort.ACCESSIONS),*copies)
            for a,b in zip(roots,copies):
                for p in (a/cohort.SNAPSHOT).iterdir():self.assertEqual(sha256(p),sha256(b/cohort.SNAPSHOT/p.name),p.name)


if __name__=='__main__':unittest.main()
