import json
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

from ds_bu_bioinformatics import module_scores
from ds_bu_bioinformatics.pca_all_cohorts import main, prepare


def cohort(root: Path, name: str, values: np.ndarray, kind=None, library_sizes=None, outcome=False, seed=0):
    rng = np.random.default_rng(seed)
    n_genes, n_samples = values.shape
    features = [f'f{i}' for i in range(n_genes)]
    samples = [f'{name}_S{j}' for j in range(n_samples)]
    folder = root / 'processed' / name
    folder.mkdir(parents=True)
    pd.DataFrame(values, index=pd.Index(features, name='feature_id'), columns=samples).to_csv(folder / 'expression.csv.gz')
    symbols = [f'A{i}' for i in range(6)] + [f'B{i}' for i in range(n_genes - 6)]
    pd.DataFrame({'feature_id': features, 'gene_id': [str(i) for i in range(n_genes)],
                  'gene_symbol': symbols}).to_csv(folder / 'genes.csv', index=False)
    clinical = pd.DataFrame({'sample_id': samples})
    if outcome:
        clinical['bcr_event'] = rng.integers(0, 2, n_samples)
        clinical['bcr_months'] = rng.uniform(5, 100, n_samples).round(1)
        clinical['bcr_eligible'] = True
    clinical.to_csv(folder / 'clinical.csv', index=False)
    expression = {'kind': kind}
    if library_sizes is not None:
        expression['library_sizes'] = dict(zip(samples, library_sizes))
    (folder / 'metadata.json').write_text(json.dumps({'expression': expression, 'primary_endpoint': 'bcr'}))
    (folder / 'validation.json').write_text(json.dumps({'structurally_valid': True}))
    return folder


def modules_file(root: Path) -> Path:
    path = root / 'modules.csv'
    rows = [{'gene': f'A{i}', 'module': 'MA'} for i in range(6)] + \
           [{'gene': g, 'module': 'TINY'} for g in ('B0', 'B1')] + [{'gene': 'ZZZ', 'module': 'ABSENT'}]
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


class PrepareTest(unittest.TestCase):
    def test_transforms_follow_declared_kind(self):
        frame = pd.DataFrame({'s1': [100.0, 900.0], 's2': [300.0, 700.0]})
        counts, how = prepare(frame, 900.0, {'expression': {'kind': 'raw_counts',
                                                             'library_sizes': {'s1': 1e6, 's2': 2e6}}}, 'x')
        np.testing.assert_allclose(counts.s2, np.log2(np.array([300, 700]) / 2e6 * 1e6 + 1))
        self.assertIn('CPM', how)
        linear, how = prepare(frame, 900.0, {'expression': {'kind': 'rpkm'}}, 'x')
        np.testing.assert_allclose(linear, np.log2(frame + 1))
        logged, how = prepare(frame / 100, 9.0, {'expression': {}}, 'x')
        np.testing.assert_allclose(logged, frame / 100)
        self.assertIn('none', how)

    def test_raw_counts_need_library_sizes_for_every_sample(self):
        frame = pd.DataFrame({'s1': [1.0], 's2': [2.0]})
        with self.assertRaises(ValueError):
            prepare(frame, 2.0, {'expression': {'kind': 'raw_counts', 'library_sizes': {'s1': 10}}}, 'x')


class EndToEndTest(unittest.TestCase):
    def test_scores_match_module_scores_and_cover_cohorts_without_outcomes(self):
        rng = np.random.default_rng(1)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            logged = rng.normal(8, 1, size=(10, 30))
            cohort(root, 'withoutcome', logged, outcome=True)
            cohort(root, 'nooutcome', rng.normal(8, 1, size=(10, 25)))
            cohort(root, 'counts', rng.integers(10, 5000, size=(10, 20)).astype(float), kind='raw_counts',
                   library_sizes=rng.integers(1_000_000, 3_000_000, size=20).tolist())
            broken = cohort(root, 'broken', rng.normal(8, 1, size=(10, 5)))
            (broken / 'genes.csv').write_text('feature_id,gene_id\nf0,0\n')
            modules = modules_file(root)
            out = root / 'out'
            with redirect_stdout(StringIO()):
                code = main(['--data-dir', str(root), '--modules', str(modules), '--output-dir', str(out)])
            self.assertEqual(code, 1)    # one cohort failed, the rest were still scored
            for name in ('withoutcome', 'nooutcome', 'counts'):
                self.assertTrue((out / f'{name}_scores.csv').is_file(), name)
            self.assertFalse((out / 'broken_scores.csv').exists())
            self.assertIn('broken', (out / 'summary.md').read_text())
            # Identical to module-scores where both apply.
            ours = pd.read_csv(out / 'withoutcome_scores.csv', index_col=0)
            theirs, _ = module_scores.score_dataset(root / 'processed' / 'withoutcome', module_scores.load_modules(modules))
            np.testing.assert_allclose(ours['MA'].to_numpy(), theirs['MA'].to_numpy(), rtol=1e-5)
            coverage = pd.read_csv(out / 'coverage.csv')
            tiny = coverage[(coverage.dataset == 'nooutcome') & (coverage.module == 'TINY')].iloc[0]
            self.assertFalse(tiny.scored)
            self.assertEqual(int(coverage[coverage.module == 'ABSENT'].genes_used.max()), 0)
            cohorts = pd.read_csv(out / 'cohorts.csv').set_index('dataset')
            self.assertEqual(cohorts.loc['nooutcome', 'eligible_outcomes'], 'none')
            self.assertEqual(cohorts.loc['withoutcome', 'eligible_outcomes'], 'bcr (30)')
            self.assertIn('CPM', cohorts.loc['counts', 'transform'])


if __name__ == '__main__':
    unittest.main()
