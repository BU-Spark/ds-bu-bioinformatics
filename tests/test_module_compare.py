"""Known expression shifts, shared feature selection, and GSVA reproducibility."""
import unittest

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.module_compare import score_expression


class ModuleCompareTest(unittest.TestCase):
    def fixture(self):
        rng = np.random.default_rng(94)
        samples = [f'S{i}' for i in range(40)]
        latent = np.linspace(-2, 2, 40)
        values = rng.normal(8, .4, (80, 40))
        values[:10] += latent
        # The same gene represented by a low-mean, unrelated probe must not
        # change either method's chosen high-mean representative.
        values = np.vstack([values, rng.normal(1, .1, 40), np.ones(40)*8])
        expression = pd.DataFrame(values, index=[f'f{i}' for i in range(82)], columns=samples)
        genes = pd.DataFrame({'feature_id': expression.index, 'gene_symbol': [f'G{i}' for i in range(80)] + ['G0', 'CONST']})
        modules = {'signal': [f'G{i}' for i in range(10)], 'tiny': ['G1', 'G2'],
                   'unmeasured': ['MISSING'], 'constant': ['CONST']}
        return expression, genes, modules, latent

    def test_known_shift_agrees_across_methods_and_small_modules_stay_missing(self):
        expression, genes, modules, latent = self.fixture()
        pca, gsva, stats, loadings, universe = score_expression(expression, genes, modules)
        self.assertEqual(pca.columns.tolist(), ['signal'])
        self.assertEqual(gsva.columns.tolist(), ['signal'])
        self.assertEqual(len(universe), 80)
        np.testing.assert_array_equal(universe.loc['G0'], expression.loc['f0'])
        self.assertGreater(pca.signal.corr(pd.Series(latent, index=pca.index)), .98)
        self.assertGreater(gsva.signal.corr(pd.Series(latent, index=gsva.index), method='spearman'), .95)
        signal = stats.set_index('module').loc['signal']
        self.assertGreater(signal.spearman_rho, .95)
        self.assertGreater(signal.pc1_correlation_with_mean_gene_z, .99)
        self.assertTrue(stats.set_index('module').loc[['tiny','unmeasured','constant'], 'spearman_rho'].isna().all())
        self.assertEqual(len(loadings), 10)

    def test_scores_are_reproducible_and_do_not_use_sample_order(self):
        expression, genes, modules, _ = self.fixture()
        a = score_expression(expression, genes, modules)
        b = score_expression(expression.iloc[:, ::-1], genes, modules)
        for i in (0, 1):
            np.testing.assert_allclose(a[i], b[i].reindex(a[i].index), atol=1e-10)

    def test_incomplete_expression_is_rejected(self):
        expression, genes, modules, _ = self.fixture()
        expression.iloc[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, 'complete finite'):
            score_expression(expression, genes, modules)


if __name__ == '__main__':
    unittest.main()
