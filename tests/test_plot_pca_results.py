import importlib.util
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np

from ds_bu_bioinformatics.pca_all_cohorts import main as score_main
from test_pca_all_cohorts import cohort, modules_file

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'plot_pca_results.py'
spec = importlib.util.spec_from_file_location('plot_pca_results', SCRIPT)
plot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plot)


class PlotPcaResultsTest(unittest.TestCase):
    def test_writes_all_four_figures_from_pca_outputs(self):
        rng = np.random.default_rng(3)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cohort(root, 'one-v1', rng.normal(8, 1, size=(10, 30)), outcome=True)
            cohort(root, 'two-v1', rng.normal(8, 1, size=(10, 25)))
            modules = modules_file(root)
            pca = root / 'pca'
            with redirect_stdout(StringIO()):
                score_main(['--data-dir', str(root), '--modules', str(modules), '--output-dir', str(pca)])
                code = plot.main(['--pca-dir', str(pca), '--data-dir', str(root), '--modules', str(modules),
                                  '--module-info', str(root / 'missing.csv'), '--cohort', 'one-v1', '--module', 'MA'])
            self.assertEqual(code, 0)
            for name in ('1_pca_explainer.png', '2_module_coverage.png', '3_pc1_variance.png'):
                self.assertGreater((pca / 'figures' / name).stat().st_size, 1000, name)

    def test_missing_cohort_skips_figures_without_failing(self):
        rng = np.random.default_rng(4)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cohort(root, 'one-v1', rng.normal(8, 1, size=(10, 20)))
            modules = modules_file(root)
            pca = root / 'pca'
            with redirect_stdout(StringIO()) as printed:
                score_main(['--data-dir', str(root), '--modules', str(modules), '--output-dir', str(pca)])
                code = plot.main(['--pca-dir', str(pca), '--data-dir', str(root), '--modules', str(modules),
                                  '--cohort', 'absent-v1'])
            self.assertEqual(code, 0)
            self.assertIn('Skipped figure 1', printed.getvalue())
            self.assertTrue((pca / 'figures' / '2_module_coverage.png').is_file())
            self.assertFalse((pca / 'figures' / '4_scores_heatmap.png').exists())

    def test_requires_pca_outputs(self):
        with tempfile.TemporaryDirectory() as temp, self.assertRaises(SystemExit):
            plot.main(['--pca-dir', temp])


if __name__ == '__main__':
    unittest.main()
