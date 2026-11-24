import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.decipher import (canonical_symbols, directions_from_matrix, load_gene_set, main,
                                           signed_score)
try:
    from test_pooled_analysis import write_cohort          # unittest discover -s tests (CI)
except ImportError:
    from tests.test_pooled_analysis import write_cohort    # python -m unittest tests.test_decipher

REPO = Path(__file__).resolve().parents[1]


class GeneSetTest(unittest.TestCase):
    def test_repository_gene_set(self):
        table = load_gene_set(REPO / 'modules' / 'decipher_erho2013.csv')
        self.assertEqual(len(table), 19)                       # 22 markers, 3 unnamed in Erho 2013 Table 2
        self.assertEqual(table.gene.nunique(), 19)
        self.assertEqual(int(table.erho_direction.notna().sum()), 15)
        self.assertEqual(set(table.module), {'decipher_erho2013'})

    def test_bad_direction_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'set.csv'
            pd.DataFrame({'gene': ['A'], 'module': ['m'], 'erho_direction': [2]}).to_csv(path, index=False)
            with self.assertRaises(ValueError):
                load_gene_set(path)

    def test_aliases(self):
        genes = pd.DataFrame({'feature_id': ['1', '2', '3'], 'gene_symbol': ['C6orf10', 'TSBP1', 'LASP1']})
        self.assertEqual(canonical_symbols(genes).gene_symbol.tolist(), ['TSBP1', 'TSBP1', 'LASP1'])
        self.assertEqual(genes.gene_symbol.tolist()[0], 'C6orf10')   # input unchanged


class ScoreTest(unittest.TestCase):
    def test_signed_score_follows_directions(self):
        rng = np.random.default_rng(0)
        latent = rng.normal(size=200)
        matrix = pd.DataFrame({'UP': latent + rng.normal(scale=0.3, size=200),
                               'DOWN': -latent + rng.normal(scale=0.3, size=200),
                               'NODIR': rng.normal(size=200)}).T
        score, used = signed_score(matrix, pd.Series({'UP': 1, 'DOWN': -1, 'NODIR': np.nan}))
        self.assertEqual(used, 2)
        self.assertGreater(np.corrcoef(score, latent)[0, 1], 0.9)
        self.assertAlmostEqual(score.std(ddof=1), 1.0)

    def test_directions_from_matrix(self):
        rng = np.random.default_rng(1)
        labels = pd.Series([1] * 100 + [0] * 100, index=[f's{i}' for i in range(200)])
        matrix = pd.DataFrame({'UP': labels + rng.normal(scale=0.5, size=200),
                               'DOWN': -labels + rng.normal(scale=0.5, size=200)}, index=labels.index).T
        table = directions_from_matrix(matrix, labels).set_index('gene')
        self.assertEqual(table.erho_direction.to_dict(), {'UP': 1, 'DOWN': -1})
        self.assertLess(table.erho_welch_p.max(), 1e-6)


class CommandTest(unittest.TestCase):
    def test_scores_cohorts_and_pools(self):
        genes = ['THBS2', 'CAMK2N1', 'IQGAP3', 'MYBPC1', 'ANO7', 'C6orf10']
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for i, name in enumerate(['tcga-prad-x', 'mortensen-x']):
                # Genes follow their Erho direction relative to the latent risk score.
                write_cohort(root, name, 150, seed=i, genes=genes, gene_effects=[1, 1, 1, -1, -1, 1])
            out = root / 'out'
            with redirect_stdout(StringIO()):
                main(['--data-dir', str(root), '--scores-dir', str(root / 'scores'), '--reference-module', 'M1',
                      '--gene-set', str(REPO / 'modules' / 'decipher_erho2013.csv'), '--output-dir', str(out)])
            per = pd.read_csv(out / 'decipher_per_cohort.csv').set_index('dataset')
            self.assertEqual(per.loc['tcga-prad-x', 'genes_pca'], 6)        # C6orf10 matched as TSBP1
            self.assertEqual(per.loc['tcga-prad-x', 'genes_signed'], 5)     # TSBP1 has no Erho direction
            self.assertTrue((per.decipher_signed_hazard_ratio_per_sd > 1).all())
            scores = pd.read_csv(out / 'tcga-prad-x_decipher_scores.csv')
            self.assertEqual(list(scores.columns), ['sample_id', 'decipher_pca', 'decipher_signed'])
            pooled = pd.read_csv(out / 'decipher_pooled.csv')
            self.assertEqual(set(pooled.terms), {'decipher_signed', 'decipher_pca', 'M1', 'decipher_signed + M1'})
            self.assertTrue((pooled.cohorts == 2).all())


if __name__ == '__main__':
    unittest.main()
