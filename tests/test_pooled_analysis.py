import json
import math
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.pooled_analysis import (fit_cox, gleason_from_grade_group, gleason_from_total,
                                                  harmonized_clinical, main, pool)


def write_cohort(root: Path, name: str, n: int, seed: int, effect: float = 0.7, genes=None, gene_effects=None):
    """Synthetic processed snapshot (+ module score file) whose hazard rises with the module score."""
    rng = np.random.default_rng(seed)
    ids = [f'{name}-{i}' for i in range(n)]
    score = rng.normal(size=n)
    months = rng.exponential(60 / np.exp(effect * score))
    event = (months < 48).astype(int)
    months = np.minimum(months, 48)
    clinical = pd.DataFrame({'sample_id': ids, 'bcr_event': event, 'bcr_months': months.round(2), 'bcr_eligible': True})
    if name.startswith('tcga-prad'):
        clinical['legacy_patient__GLEASON_PATTERN_PRIMARY'] = rng.choice([3, 4], n)
        clinical['legacy_patient__GLEASON_PATTERN_SECONDARY'] = rng.choice([3, 4, 5], n)
        clinical['path_t_group'] = rng.choice(['T2', 'T3'], n)
    elif name.startswith('mortensen'):
        clinical['geo__gleason_grade'] = rng.choice([6, 7, 8, 9], n)
        clinical['geo__t_stage'] = rng.choice(['pT2c', 'pT3a'], n)
    elif name.startswith('sboner'):
        clinical['gleason_score'] = rng.choice([6, 7, 8], n)
    clinical.loc[0, 'bcr_eligible'] = False
    folder = root / 'processed' / name
    folder.mkdir(parents=True)
    clinical.to_csv(folder / 'clinical.csv', index=False)
    (folder / 'metadata.json').write_text(json.dumps({'primary_endpoint': 'bcr'}))
    genes = genes or [f'G{i}' for i in range(6)]
    gene_effects = gene_effects if gene_effects is not None else [1.0] * len(genes)
    expression = pd.DataFrame({g: 8 + e * score + rng.normal(scale=0.5, size=n) for g, e in zip(genes, gene_effects)},
                              index=ids).T
    expression.index = [f'f{i}' for i in range(len(genes))]
    expression.index.name = 'feature_id'
    expression.to_csv(folder / 'expression.csv.gz')
    pd.DataFrame({'feature_id': expression.index, 'gene_id': range(len(genes)), 'gene_symbol': genes}).to_csv(
        folder / 'genes.csv', index=False)
    scores_dir = root / 'scores'
    scores_dir.mkdir(exist_ok=True)
    z = (score - score.mean()) / score.std(ddof=1)
    pd.DataFrame({'sample_id': ids, 'M1': z, 'M2': rng.normal(size=n)}).to_csv(scores_dir / f'{name}_scores.csv', index=False)
    return folder


class ClinicalTest(unittest.TestCase):
    def test_gleason_levels(self):
        self.assertEqual(gleason_from_grade_group(pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])).tolist(), [0, 1, 1, 2, 2])
        levels = gleason_from_total(pd.Series(['5', '6', '7', '8', '10', 'unknown']))
        self.assertEqual(levels.tolist()[:5], [0, 0, 1, 2, 2])
        self.assertTrue(math.isnan(levels.iloc[5]))
        with self.assertRaises(ValueError):
            gleason_from_total(pd.Series(['12']))

    def test_unmapped_cohort_and_missing_stage(self):
        clinical = pd.DataFrame({'gleason_score': ['6', '8']})
        self.assertTrue(harmonized_clinical(clinical, 'belfast-v1').isna().all().all())
        sboner = harmonized_clinical(clinical, 'sboner-v2')
        self.assertEqual(sboner.gleason_cat.tolist(), [0, 2])
        self.assertTrue(sboner.pt3_plus.isna().all())


class PooledModelTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        for i, name in enumerate(['tcga-prad-x', 'mortensen-x', 'sboner-x']):
            write_cohort(self.root, name, 150, seed=i)

    def tearDown(self):
        self.tmp.cleanup()

    def test_pool_keeps_eligible_patients_and_recovers_effect(self):
        pooled = pool(self.root, self.root / 'scores', ['tcga-prad-x', 'mortensen-x', 'sboner-x'])
        self.assertEqual(len(pooled), 3 * 149)
        self.assertEqual(set(pooled.study), {'tcga-prad-x', 'mortensen-x', 'sboner-x'})
        stratified = fit_cox(pooled, ['M1'], [], 'study')
        self.assertGreater(stratified['M1_hr'], 1.5)
        self.assertLess(stratified['M1_p'], 1e-6)
        dummies = fit_cox(pooled, ['M1'], ['study'], None)
        self.assertAlmostEqual(dummies['M1_hr'], stratified['M1_hr'], delta=0.3)
        # Sboner has no T stage, so the stage-adjusted model drops that cohort.
        adjusted = fit_cox(pooled, ['M1'], ['gleason_cat', 'pt3_plus'], 'study')
        self.assertEqual(adjusted['cohorts'], 2)
        self.assertEqual(adjusted['n'], 2 * 149)

    def test_cli_writes_tables_and_orientation(self):
        modules = self.root / 'modules.csv'
        pd.DataFrame({'gene': [f'G{i}' for i in range(6)], 'module': 'M1'}).to_csv(modules, index=False)
        out = self.root / 'out'
        with redirect_stdout(StringIO()):
            main(['--data-dir', str(self.root), '--scores-dir', str(self.root / 'scores'), '--dataset', 'tcga-prad-x',
                  '--dataset', 'mortensen-x', '--modules', str(modules), '--output-dir', str(out)])
        table = pd.read_csv(out / 'pooled_cox.csv')
        self.assertEqual(set(table.model), {'stratified', 'study_covariate', 'strat_adj_gleason', 'strat_adj_gleason_pt3'})
        m1 = table[(table.module == 'M1') & (table.model == 'stratified')].iloc[0]
        self.assertLess(m1.q, 0.05)
        self.assertTrue((out / 'pooled_cox_wide.csv').exists())
        orientation = pd.read_csv(out / 'pc1_orientation.csv')
        self.assertTrue((orientation.corr_with_mean_expression > 0).all())


if __name__ == '__main__':
    unittest.main()
