"""Pooled patient-level Cox across cohorts with a study term, plus a PC1 orientation check.

Complements `module-analysis` (per-cohort Cox, then a random-effects meta-analysis) with the model
requested at the 11/16 client meeting: one Cox model over the patients of several cohorts with a
term for study, so differences between datasets are accounted for.

  * Scores come from `module-scores` (<snapshot>_scores.csv) and are already standardized within
    each cohort (mean 0, SD 1), so hazard ratios stay "per SD of score".
  * Each cohort keeps its own endpoint (PFI, BCR, DFS, cancer death). The primary model is
    stratified by study, so every cohort has its own baseline hazard; a model with study as dummy
    covariates is reported as well.
  * Adjusted models add Gleason (three levels: <=6, 7, 8+, so cohorts that record only the total
    score can be included) and pathological T3+. Patients missing a covariate are dropped from that
    model, which removes whole cohorts that do not record it (for example Sboner has no stage).
  * Benjamini-Hochberg q-values are computed across modules within each model.

The orientation check re-derives each module score and reports its correlation with the mean
standardized expression of the module genes. The scoring rule flips PC1 so this correlation is
positive (higher score = higher module expression); a low value means the module's genes are not
co-expressed in that cohort, so PC1 summarizes them poorly.

Outputs (in --output-dir): pooled_cox.csv (long), pooled_cox_wide.csv, pooled_cohorts.csv and,
with --modules, pc1_orientation.csv.
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.module_analysis import bh, grade_from_patterns, grade_from_text, number, pt3_plus
from ds_bu_bioinformatics.module_scores import endpoint, load_dataset, load_modules, module_matrix, pca_score

# Cohorts in the five-cohort summary presented on 11/16; Taylor and Long are provisional.
RELIABLE = ('tcga-prad', 'stockholm', 'mortensen', 'sboner', 'cambridge')
LOW_COHERENCE = 0.5


def gleason_from_grade_group(grade_group: pd.Series) -> pd.Series:
    """0: Gleason <=6 (GG1), 1: Gleason 7 (GG2-3), 2: Gleason 8-10 (GG4-5)."""
    return grade_group.map({1.0: 0.0, 2.0: 1.0, 3.0: 1.0, 4.0: 2.0, 5.0: 2.0})


def gleason_from_total(total: pd.Series) -> pd.Series:
    """Same three levels from a total Gleason score, for cohorts without primary/secondary patterns."""
    t = number(total)
    if ((t < 2) | (t > 10)).any():
        raise ValueError(f'Gleason score outside 2-10: {sorted(t[(t < 2) | (t > 10)].unique())}')
    return pd.Series(np.select([t <= 6, t == 7], [0.0, 1.0], 2.0), index=total.index).where(t.notna())


def _none(clinical: pd.DataFrame) -> pd.Series:
    return pd.Series(np.nan, index=clinical.index)


# Snapshot name prefix -> (gleason, pt3_plus) from that cohort's processed clinical.csv.
CLINICAL = {
    'tcga-prad': (lambda c: gleason_from_grade_group(grade_from_patterns(c['legacy_patient__GLEASON_PATTERN_PRIMARY'],
                                                                         c['legacy_patient__GLEASON_PATTERN_SECONDARY'])),
                  lambda c: pt3_plus(c['path_t_group'])),
    'cambridge': (lambda c: gleason_from_grade_group(grade_from_text(c['geo__tumour_gleason'])),
                  lambda c: pt3_plus(c['geo__pathology_stage'])),
    'stockholm': (lambda c: gleason_from_grade_group(grade_from_text(c['geo__tumour_gleason'])),
                  lambda c: pt3_plus(c['geo__pathology_stage'])),
    'mortensen': (lambda c: gleason_from_total(c['geo__gleason_grade']), lambda c: pt3_plus(c['geo__t_stage'])),
    # Watchful-waiting TURP cohort: GEO records a cohort-level stage only, so no per-patient T stage.
    'sboner': (lambda c: gleason_from_total(c['gleason_score']), _none),
}


def harmonized_clinical(clinical: pd.DataFrame, dataset: str) -> pd.DataFrame:
    """gleason_cat and pt3_plus for one cohort; all missing when the cohort has no mapping."""
    spec = next((v for k, v in CLINICAL.items() if dataset.startswith(k)), None)
    if spec is None:
        return pd.DataFrame({'gleason_cat': np.nan, 'pt3_plus': np.nan}, index=clinical.index)
    gleason, stage = spec
    return pd.DataFrame({'gleason_cat': gleason(clinical).astype(float), 'pt3_plus': stage(clinical).astype(float)},
                        index=clinical.index)


def cohort_frame(processed: Path, scores: pd.DataFrame) -> pd.DataFrame:
    """Eligible patients of one cohort: study, endpoint, months, event, clinical covariates and score columns.

    `scores` is indexed by sample_id in clinical.csv order (as module-scores writes it)."""
    name = processed.name
    clinical = pd.read_csv(processed / 'clinical.csv', dtype={'sample_id': str}, low_memory=False)
    metadata = json.loads((processed / 'metadata.json').read_text())
    if scores.index.tolist() != clinical.sample_id.tolist():
        raise ValueError(f'{name}: score samples do not match clinical.csv; rerun module-scores')
    endpoint_name, event, months, eligible = endpoint(clinical, metadata, name)
    covariates = harmonized_clinical(clinical, name)
    frame = pd.DataFrame({'study': name, 'endpoint': endpoint_name, 'months': months.to_numpy(float),
                          'event': event.to_numpy(float), 'gleason_cat': covariates.gleason_cat.to_numpy(),
                          'pt3_plus': covariates.pt3_plus.to_numpy()}, index=clinical.sample_id)
    return frame.join(scores)[eligible]


def pool(data_dir: Path, scores_dir: Path, datasets) -> pd.DataFrame:
    frames = []
    for name in datasets:
        scores = pd.read_csv(scores_dir / f'{name}_scores.csv', dtype={'sample_id': str}).set_index('sample_id')
        frames.append(cohort_frame(data_dir / 'processed' / name, scores))
    return pd.concat(frames)


# model -> (adjustment covariates, strata). 'study' as a covariate is expanded to dummies.
MODELS = {
    'stratified': ([], 'study'),
    'study_covariate': (['study'], None),
    'strat_adj_gleason': (['gleason_cat'], 'study'),
    'strat_adj_gleason_pt3': (['gleason_cat', 'pt3_plus'], 'study'),
}


def design(frame: pd.DataFrame, terms: list[str], covariates: list[str], strata: str | None) -> pd.DataFrame:
    """Complete-case design matrix: terms, dummy-coded covariates, months, event (+ strata column)."""
    needed = terms + [c for c in covariates if c != 'study'] + ['months', 'event']
    frame = frame.dropna(subset=needed)
    parts = [frame[terms]]
    for c in covariates:
        if c in ('study', 'gleason_cat'):
            parts.append(pd.get_dummies(frame[c], prefix=c, drop_first=True, dtype=float))
        else:
            parts.append(frame[[c]])
    x = pd.concat(parts, axis=1)
    x = x.loc[:, [c for c in x if c in terms or x[c].nunique() > 1]]
    x[['months', 'event']] = frame[['months', 'event']]
    if strata:
        x[strata] = frame[strata]
    return x


def fit_cox(frame: pd.DataFrame, terms: list[str], covariates=(), strata: str | None = None) -> dict:
    """Cox fit; returns hazard ratio, 95% CI and p for each term plus n, events, cohorts and C-index."""
    from lifelines import CoxPHFitter
    x = design(frame, list(terms), list(covariates), strata)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fit = CoxPHFitter().fit(x, duration_col='months', event_col='event', strata=[strata] if strata else None)
    out = {'n': len(x), 'events': int(x.event.sum()), 'cohorts': int(frame.loc[x.index, 'study'].nunique()),
           'c_index': float(fit.concordance_index_)}
    for t in terms:
        row = fit.summary.loc[t]
        out.update({f'{t}_hr': float(row['exp(coef)']), f'{t}_ci_low': float(row['exp(coef) lower 95%']),
                    f'{t}_ci_high': float(row['exp(coef) upper 95%']), f'{t}_p': float(row['p'])})
    return out


def pooled_table(pooled: pd.DataFrame, score_columns, models=None) -> pd.DataFrame:
    """One row per score and model, with BH q-values across scores within each model."""
    rows = []
    for column in score_columns:
        for model, (covariates, strata) in (models or MODELS).items():
            record = {'module': column, 'model': model}
            try:
                result = fit_cox(pooled, [column], covariates, strata)
                record.update({k.replace(f'{column}_', ''): v for k, v in result.items()})
            except Exception as error:   # convergence failures are reported, not fatal
                record['note'] = f'Cox fit failed: {error}'
            rows.append(record)
    table = pd.DataFrame(rows)
    table['q'] = table.groupby('model').p.transform(lambda p: pd.Series(bh(p), index=p.index))
    return table


def orientation_check(data_dir: Path, scores_dir: Path, modules: dict, datasets) -> pd.DataFrame:
    rows = []
    for name in datasets:
        expression, genes, *_ = load_dataset(data_dir / 'processed' / name)
        saved = pd.read_csv(scores_dir / f'{name}_scores.csv', dtype={'sample_id': str}).set_index('sample_id')
        for module, symbols in modules.items():
            if module not in saved:
                continue
            matrix = module_matrix(expression, genes, symbols)
            z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=1), axis=0)
            score, variance, loadings = pca_score(matrix)
            stored = saved[module].reindex(score.index)
            rows.append({'dataset': name, 'module': module, 'genes_used': len(matrix),
                         'corr_with_mean_expression': float(np.corrcoef(stored, z.mean(axis=0))[0, 1]),
                         'matches_stored_score': bool(np.allclose(stored, score, atol=1e-6)),
                         'share_positive_loadings': float((loadings > 0).mean()), 'pc1_variance': variance})
    table = pd.DataFrame(rows)
    table['low_coherence'] = table.corr_with_mean_expression < LOW_COHERENCE
    return table


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--data-dir', type=Path, default=Path('data'), help='Folder containing processed/<snapshot>')
    parser.add_argument('--scores-dir', type=Path, help='module-scores output; default <data-dir>/analysis/module_scores')
    parser.add_argument('--dataset', action='append', default=[],
                        help='Snapshot to pool (repeatable). Default: the five reliable cohorts found in --scores-dir')
    parser.add_argument('--modules', type=Path, help='gene,module CSV; enables the PC1 orientation check')
    parser.add_argument('--module-info', type=Path, help='module,category,annotation CSV for readable labels')
    parser.add_argument('--output-dir', type=Path, help='Default: <scores-dir>/pooled')
    args = parser.parse_args(argv)
    scores_dir = args.scores_dir or args.data_dir / 'analysis' / 'module_scores'
    output = args.output_dir or scores_dir / 'pooled'
    output.mkdir(parents=True, exist_ok=True)
    datasets = args.dataset or sorted(p.name[:-len('_scores.csv')] for p in scores_dir.glob('*_scores.csv')
                                      if p.name.startswith(RELIABLE))
    if not datasets:
        parser.error(f'no score files for {RELIABLE} in {scores_dir}')

    pooled = pool(args.data_dir, scores_dir, datasets)
    score_columns = [c for c in pooled.columns
                     if c not in {'study', 'endpoint', 'months', 'event', 'gleason_cat', 'pt3_plus'}]
    cohorts = pooled.groupby('study').agg(endpoint=('endpoint', 'first'), patients=('event', 'size'),
                                          events=('event', 'sum'),
                                          gleason_missing=('gleason_cat', lambda s: int(s.isna().sum())),
                                          pt3_missing=('pt3_plus', lambda s: int(s.isna().sum())))
    cohorts.to_csv(output / 'pooled_cohorts.csv')
    table = pooled_table(pooled, score_columns)
    if args.module_info:
        info = pd.read_csv(args.module_info)[['module', 'category', 'annotation']]
        table = table.merge(info, on='module', how='left')
    table.to_csv(output / 'pooled_cox.csv', index=False)
    wide = table.pivot(index='module', columns='model', values=['hr', 'ci_low', 'ci_high', 'p', 'q'])
    wide.columns = [f'{model}_{value}' for value, model in wide.columns]
    wide.sort_values('stratified_p').to_csv(output / 'pooled_cox_wide.csv')

    print(cohorts.to_string(), '\n')
    for model in MODELS:
        m = table[table.model == model]
        sig = m[m.q < 0.05]
        print(f'{model:24s} q < 0.05: {len(sig):3d} of {int(m.p.notna().sum())} '
              f'(higher risk {int((sig.hr > 1).sum())}, lower risk {int((sig.hr < 1).sum())}); '
              f'patients {int(m.n.max())}, cohorts {int(m.cohorts.max())}')
    if args.modules:
        orientation = orientation_check(args.data_dir, scores_dir, load_modules(args.modules), datasets)
        orientation.to_csv(output / 'pc1_orientation.csv', index=False)
        print(f'\nPC1 orientation: {len(orientation)} scores, stored scores reproduced: '
              f'{orientation.matches_stored_score.mean():.0%}, positive correlation with mean expression: '
              f'{(orientation.corr_with_mean_expression > 0).mean():.0%}, '
              f'low coherence (< {LOW_COHERENCE}): {int(orientation.low_coherence.sum())}')
    print(f'\nWrote {output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
