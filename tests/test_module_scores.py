import json
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.module_scores import load_modules, main, module_matrix, pca_score, score_dataset, survival
from ds_bu_bioinformatics.analysis_policy import survival_blocker

MODULES = {'prolif': [f'P{i}' for i in range(10)], 'tiny': ['P0', 'P1', 'ABSENT']}


def write_dataset(root: Path, name='demo-v1', scale=1.0, endpoint_style='shared'):
    """Synthetic dataset: module genes follow a latent factor that also shortens time to event."""
    rng = np.random.default_rng(1)
    n = 80
    samples = [f'S{i:02d}' for i in range(n)]
    factor = rng.normal(size=n)
    rows, symbols = [], []
    for i in range(10):
        rows.append(8 + factor + rng.normal(scale=0.5, size=n))
        symbols.append(f'P{i}')
    rows.append(3 + rng.normal(scale=0.5, size=n))   # second, lower-mean probe for P0: must not be used
    symbols.append('P0')
    for i in range(5):
        rows.append(6 + rng.normal(size=n))
        symbols.append(f'N{i}')
    feature_ids = [f'f{i}' for i in range(len(rows))]
    values = np.array(rows)
    if scale != 1.0:
        values = 2 ** values * scale        # unlogged scale, like Taylor
    processed = root / 'processed' / name
    processed.mkdir(parents=True)
    pd.DataFrame(values, index=pd.Index(feature_ids, name='feature_id'), columns=samples).to_csv(
        processed / 'expression.csv.gz')
    pd.DataFrame({'feature_id': feature_ids, 'gene_id': [str(i) for i in range(len(rows))],
                  'gene_symbol': symbols}).to_csv(processed / 'genes.csv', index=False)
    months = np.round(rng.exponential(scale=30 * np.exp(-factor)), 2) + 0.1
    event = (rng.uniform(size=n) < 0.7).astype(int)
    clinical = pd.DataFrame({'sample_id': samples})
    if endpoint_style == 'shared':
        clinical = clinical.assign(bcr_event=event, bcr_months=months, bcr_eligible=True,
                                   bcr_eligibility_reason='eligible', bcr_time_unit='months')
    else:   # Sboner-style: one follow-up time and eligibility flag for every endpoint
        clinical = clinical.assign(pcsm_event=event, os_event=event, followup_months=months, survival_eligible=True)
    clinical.to_csv(processed / 'clinical.csv', index=False)
    (processed / 'metadata.json').write_text(json.dumps({'snapshot': name}))
    return processed, factor


class ModuleScoreTest(unittest.TestCase):
    def test_no_comparable_survival_pairs_does_not_abort_scoring(self):
        rng = np.random.default_rng(17)
        result = survival(pd.Series(rng.normal(size=25)),
                          pd.Series([10.] * 20 + [1.] * 5), pd.Series([1] * 20 + [0] * 5))
        self.assertIn('hazard_ratio_per_sd', result)
        self.assertTrue(np.isnan(result['concordance']))
        self.assertIn('Concordance unavailable', result['note'])

    def test_sboner_extreme_groups_block_cox_but_keep_pca_scores(self):
        self.assertIn('Outcome-selected', survival_blocker({'geo_series': 'GSE16560'}))
        processed, _ = write_dataset(self.root, name='sboner-v1', endpoint_style='sboner')
        (processed / 'metadata.json').write_text(json.dumps({'geo_series': 'GSE16560'}))
        scores, summary = score_dataset(processed, MODULES)
        self.assertEqual(len(scores), 80)
        self.assertTrue(summary.hazard_ratio_per_sd.isna().all())
        self.assertIn('Outcome-selected', summary.note.iloc[0])

    def test_status_only_cohort_still_writes_pca_scores_without_survival(self):
        processed, factor = write_dataset(self.root)
        clinical = pd.read_csv(processed / 'clinical.csv')
        clinical['bcr_months'] = np.nan
        clinical['bcr_eligible'] = False
        clinical.to_csv(processed / 'clinical.csv', index=False)
        scores, summary = score_dataset(processed, MODULES)
        self.assertGreater(np.corrcoef(scores.prolif, factor)[0, 1], 0.9)
        self.assertTrue(summary.hazard_ratio_per_sd.isna().all())
        self.assertEqual(summary.eligible.tolist(), [0, 0])
        self.assertIn('PCA scored; survival not fitted', summary.note.iloc[0])
        # A whole run with no fitted models must still write a readable report/plot.
        modules = self.root / 'modules.csv'
        modules.write_text('gene,module\n' + ''.join(f'{g},{m}\n' for m, gs in MODULES.items() for g in gs))
        with redirect_stdout(StringIO()):
            main(['--data-dir', str(self.root), '--modules', str(modules)])
        out = self.root / 'analysis' / 'module_scores'
        self.assertTrue((out / 'hazard_ratios.png').is_file())
        self.assertIn('PCA scored; survival not fitted', (out / 'summary.md').read_text())
        self.assertEqual(pd.read_csv(out / 'demo-v1_scores.csv').shape, (80, 2))

    def test_all_censored_cohort_does_not_attempt_a_cox_fit(self):
        processed, _ = write_dataset(self.root)
        clinical = pd.read_csv(processed / 'clinical.csv')
        clinical['bcr_event'] = 0
        clinical.to_csv(processed / 'clinical.csv', index=False)
        scores, summary = score_dataset(processed, MODULES)
        self.assertEqual(len(scores), 80)
        self.assertTrue(summary.hazard_ratio_per_sd.isna().all())

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_one_feature_per_gene_highest_mean(self):
        processed, _ = write_dataset(self.root)
        expression = pd.read_csv(processed / 'expression.csv.gz', index_col=0)
        genes = pd.read_csv(processed / 'genes.csv', dtype=str)
        matrix = module_matrix(expression, genes, MODULES['prolif'])
        self.assertEqual(len(matrix), 10)
        np.testing.assert_array_equal(matrix.loc['P0'].to_numpy(), expression.loc['f0'].to_numpy())

    def test_pc1_sign_follows_module_expression(self):
        processed, factor = write_dataset(self.root)
        expression = pd.read_csv(processed / 'expression.csv.gz', index_col=0)
        genes = pd.read_csv(processed / 'genes.csv', dtype=str)
        score, variance, loadings = pca_score(module_matrix(expression, genes, MODULES['prolif']))
        self.assertGreater(np.corrcoef(score, factor)[0, 1], 0.9)
        self.assertGreater(variance, 0.5)
        self.assertTrue((loadings > 0).all())
        self.assertAlmostEqual(score.std(ddof=1), 1.0)
        self.assertAlmostEqual(score.mean(), 0.0)
        # Flipping the data must not flip the score's meaning.
        flipped, _, _ = pca_score(-module_matrix(expression, genes, MODULES['prolif']))
        self.assertLess(np.corrcoef(flipped, factor)[0, 1], -0.9)

    def test_score_dataset_reports_hazard_and_skips_small_modules(self):
        processed, _ = write_dataset(self.root)
        scores, summary = score_dataset(processed, MODULES)
        self.assertEqual(scores.columns.tolist(), ['prolif'])
        row = summary.set_index('module').loc['prolif']
        self.assertEqual((row.endpoint, row.eligible, row.genes_used), ('bcr', 80, 10))
        self.assertGreater(row.hazard_ratio_per_sd, 1.5)
        self.assertLess(row.p_value, 0.001)
        self.assertGreater(row.concordance, 0.6)
        self.assertIn('not scored', summary.set_index('module').loc['tiny', 'note'])

    def test_unlogged_expression_is_log_transformed(self):
        processed, factor = write_dataset(self.root, scale=10.0)
        scores, summary = score_dataset(processed, MODULES)
        self.assertTrue(summary.expression_transform.str.startswith('log2').all())
        self.assertGreater(np.corrcoef(scores.prolif, factor)[0, 1], 0.9)

    def test_concordance_uses_fitted_risk_direction_for_protective_modules(self):
        processed, factor = write_dataset(self.root)
        clinical = pd.read_csv(processed / 'clinical.csv')
        score = pd.Series(factor)
        harmful = survival(score, clinical.bcr_months, clinical.bcr_event)
        protective = survival(-score, clinical.bcr_months, clinical.bcr_event)
        self.assertGreater(harmful['hazard_ratio_per_sd'], 1.0)
        self.assertLess(protective['hazard_ratio_per_sd'], 1.0)
        self.assertGreater(protective['concordance'], 0.6)
        self.assertAlmostEqual(harmful['concordance'], protective['concordance'])

    def test_shared_followup_fallback_for_sboner(self):
        processed, _ = write_dataset(self.root, name='sboner-v1', endpoint_style='sboner')
        _, summary = score_dataset(processed, MODULES)
        self.assertEqual(summary.endpoint.iloc[0], 'pcsm')
        processed, _ = write_dataset(self.root, name='other-v1', endpoint_style='sboner')
        with self.assertRaisesRegex(ValueError, 'primary endpoint'):
            score_dataset(processed, MODULES)

    def test_cli_writes_summary_scores_and_plot(self):
        write_dataset(self.root)
        modules = self.root / 'modules.csv'
        modules.write_text('gene,module\n' + ''.join(f'{g},{m}\n' for m, gs in MODULES.items() for g in gs))
        self.assertEqual(sorted(load_modules(modules)), ['prolif', 'tiny'])
        with redirect_stdout(StringIO()) as output:
            main(['--data-dir', str(self.root), '--modules', str(modules)])
        out = self.root / 'analysis' / 'module_scores'
        for name in ('summary.csv', 'summary.md', 'hazard_ratios.png', 'demo-v1_scores.csv'):
            self.assertTrue((out / name).is_file(), name)
        self.assertIn('| demo-v1 | BCR (80,', output.getvalue())


if __name__ == '__main__':
    unittest.main()
