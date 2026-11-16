import json
import math
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.module_analysis import (bh, clinical_covariates, gleason_8_plus, grade_from_patterns,
                                                  grade_from_text, main, n1, pt3_plus, random_effects)


class ParserTest(unittest.TestCase):
    def test_grade_group_from_text_and_patterns(self):
        text = pd.Series(['6=3+3', '7=3+4', '7=4+3', '8=4+4', '9=5+4', 'unknown'])
        self.assertEqual(grade_from_text(text).tolist()[:5], [1, 2, 3, 4, 5])
        self.assertTrue(math.isnan(grade_from_text(text).iloc[5]))
        primary, secondary = pd.Series(['3', '4', '4.0', None]), pd.Series(['4', '3', '5', '3'])
        result = grade_from_patterns(primary, secondary)
        self.assertEqual(result.tolist()[:3], [2, 3, 5])
        self.assertTrue(math.isnan(result.iloc[3]))

    def test_bad_values_fail_instead_of_becoming_missing(self):
        with self.assertRaises(ValueError):
            grade_from_text(pd.Series(['7=3+3']))       # patterns do not add up
        with self.assertRaises(ValueError):
            pt3_plus(pd.Series(['stage two']))
        with self.assertRaises(ValueError):
            n1(pd.Series(['N2']))

    def test_t_and_n_stage_parsing(self):
        stages = pt3_plus(pd.Series(['pT3a N0Mx', 'T2bMxNx', 'T4', 'T2C', 'T0MxNx', 'TxM0N0', 'UNKNOWN']))
        self.assertEqual(stages.tolist()[:5], [1, 0, 1, 0, 0])
        self.assertTrue(stages.iloc[5:].isna().all())
        self.assertEqual(n1(pd.Series(['N0', 'N1'])).tolist(), [0, 1])
        self.assertTrue(math.isnan(n1(pd.Series(['NX'])).iloc[0]))

    def test_belfast_total_gleason_and_clinical_stage(self):
        self.assertEqual(gleason_8_plus(pd.Series(['6', '7', '8', '10'])).tolist(), [0, 0, 1, 1])
        self.assertTrue(math.isnan(gleason_8_plus(pd.Series(['unknown'])).iloc[0]))
        with self.assertRaises(ValueError):
            gleason_8_plus(pd.Series(['12']))
        clinical = pd.DataFrame({'geo__patient_age_years': ['70', '65'], 'geo__gleason_grade': ['9', '7'],
                                 'geo__t_stage': ['T3', 'Unknown']})
        covariates = clinical_covariates(clinical, 'belfast-v1')
        self.assertEqual(list(covariates), ['age', 'gleason_8_plus', 'ct3_plus'])
        self.assertEqual(covariates.iloc[0].tolist(), [70, 1, 1])
        self.assertTrue(math.isnan(covariates.ct3_plus.iloc[1]))

    def test_unmapped_cohort_has_no_covariates(self):
        self.assertTrue(clinical_covariates(pd.DataFrame({'x': [1]}), 'sboner-v1').empty)


class StatisticsTest(unittest.TestCase):
    def test_bh_matches_hand_calculation_and_keeps_missing(self):
        q = bh([0.01, 0.04, 0.03, np.nan])
        np.testing.assert_allclose(q[:3], [0.03, 0.04, 0.04])
        self.assertTrue(np.isnan(q[3]))

    def test_identical_studies_pool_without_heterogeneity(self):
        result = random_effects(np.log([1.5, 1.5, 1.5]), [0.2, 0.2, 0.2])
        self.assertAlmostEqual(result['hr'], 1.5)
        self.assertEqual(result['i2'], 0.0)
        self.assertEqual(result['same_direction'], 3)
        self.assertAlmostEqual(math.log(result['ci_high']) - math.log(1.5), 1.959963984540054 * 0.2 / math.sqrt(3))

    def test_conflicting_studies_show_heterogeneity(self):
        result = random_effects(np.log([2.0, 0.5]), [0.1, 0.1])
        self.assertAlmostEqual(result['hr'], 1.0)
        self.assertGreater(result['i2'], 0.9)
        self.assertGreater(result['tau2'], 0)
        self.assertEqual(random_effects([0.1], [0.1]), {'k': 1})


def synthetic_cohort(root: Path, name: str, seed: int, n: int = 120):
    """A processed snapshot plus module-scores output in which module A tracks risk and B is noise."""
    rng = np.random.default_rng(seed)
    samples = [f'{name}_S{i:03}' for i in range(n)]
    score_a, score_b = rng.normal(size=n), rng.normal(size=n)
    grade = rng.integers(1, 6, size=n)
    hazard = np.exp(0.7 * score_a + 0.2 * grade)
    time = rng.exponential(60 / hazard)
    censor = rng.uniform(20, 120, size=n)
    gleason = {1: '6=3+3', 2: '7=3+4', 3: '7=4+3', 4: '8=4+4', 5: '9=4+5'}
    clinical = pd.DataFrame({'sample_id': samples, 'bcr_event': (time <= censor).astype(int),
                             'bcr_months': np.minimum(time, censor).round(2), 'bcr_eligible': True,
                             'geo__tumour_gleason': [gleason[g] for g in grade],
                             'geo__pathology_stage': rng.choice(['pT2c N0Mx', 'pT3a N0Mx'], size=n),
                             'geo__age_at_diag': rng.integers(50, 75, size=n).astype(float)})
    clinical.loc[0, 'geo__tumour_gleason'] = 'unknown'          # one incomplete case
    processed = root / 'processed' / name
    processed.mkdir(parents=True)
    clinical.to_csv(processed / 'clinical.csv', index=False)
    (processed / 'metadata.json').write_text(json.dumps({'primary_endpoint': 'bcr'}))
    scores_dir = root / 'analysis' / 'module_scores'
    scores_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({'sample_id': samples, 'A': score_a, 'B': score_b}).to_csv(scores_dir / f'{name}_scores.csv', index=False)
    rows = []
    for module, hr, p in (('A', 2.0, 1e-6), ('B', 1.0, 0.9)):
        rows.append({'dataset': name, 'module': module, 'endpoint': 'bcr', 'hazard_ratio_per_sd': hr,
                     'hr_ci_low': hr / 1.4, 'hr_ci_high': hr * 1.4, 'p_value': p})
    rows.append({'dataset': name, 'module': 'tiny', 'endpoint': 'bcr', 'note': 'not scored'})
    return rows


class EndToEndTest(unittest.TestCase):
    def test_cli_writes_adjusted_and_meta_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rows = synthetic_cohort(root, 'cambridge-v9', 1) + synthetic_cohort(root, 'stockholm-v9', 2)
            scores_dir = root / 'analysis' / 'module_scores'
            pd.DataFrame(rows).to_csv(scores_dir / 'summary.csv', index=False)
            info = root / 'info.csv'
            pd.DataFrame({'module': ['A', 'B'], 'category': ['x', 'y'], 'annotation': ['Risky', '']}).to_csv(info, index=False)
            with redirect_stdout(StringIO()) as printed:
                self.assertEqual(main(['--data-dir', str(root), '--module-info', str(info), '--top', '5']), 0)
            out = scores_dir / 'analysis'
            for name in ('per_dataset.csv', 'meta.csv', 'clinical_models.csv', 'analysis.md', 'heatmap.png', 'forest_top.png'):
                self.assertTrue((out / name).is_file(), name)
            table = pd.read_csv(out / 'per_dataset.csv')
            a = table[(table.module == 'A') & (table.dataset == 'cambridge-v9')].iloc[0]
            self.assertGreater(a.adj_hazard_ratio_per_sd, 1.3)          # simulated effect survives adjustment
            self.assertLess(a.added_value_lr_p, 0.05)
            self.assertTrue(table[table.module == 'tiny'].adj_hazard_ratio_per_sd.isna().all())
            models = pd.read_csv(out / 'clinical_models.csv').set_index('dataset')
            self.assertEqual(models.loc['cambridge-v9', 'covariates_used'], 'age, grade_group, pt3_plus')
            self.assertEqual(models.loc['cambridge-v9', 'complete_cases'], 119)
            meta = pd.read_csv(out / 'meta.csv')
            self.assertEqual(meta.iloc[0].module, 'A')
            self.assertEqual(meta.iloc[0].label, 'A · Risky')
            self.assertEqual(int(meta.iloc[0].k), 2)
            self.assertIn('A · Risky', printed.getvalue())

    def test_misaligned_scores_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rows = synthetic_cohort(root, 'stockholm-v9', 3)
            scores_dir = root / 'analysis' / 'module_scores'
            pd.DataFrame(rows).to_csv(scores_dir / 'summary.csv', index=False)
            path = scores_dir / 'stockholm-v9_scores.csv'
            pd.read_csv(path).iloc[::-1].to_csv(path, index=False)
            with self.assertRaises(ValueError):
                main(['--data-dir', str(root)])


if __name__ == '__main__':
    unittest.main()
