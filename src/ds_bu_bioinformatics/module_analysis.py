"""Statistics on top of module-scores output: FDR, clinical-adjusted Cox and cross-cohort meta-analysis.

Reads what `module-scores` wrote (summary.csv and <snapshot>_scores.csv) plus each snapshot's
processed clinical.csv. Nothing upstream is changed. For every dataset and module:

  1. Unadjusted results are taken from module-scores' summary.csv, and Benjamini-Hochberg q-values
     are added within each dataset.
  2. A Cox model of the dataset's own endpoint on score + clinical covariates (complete cases among
     eligible records). Covariates are mapped per cohort in COVARIATES because each source names and
     codes them differently. The added value of the score over a clinical-only model on the same
     patients is reported as a likelihood-ratio test and a change in C-index.
  3. A DerSimonian-Laird random-effects meta-analysis across datasets of the log hazard ratio per SD,
     separately for unadjusted and adjusted estimates, with I² heterogeneity and BH q-values across
     modules. Endpoints differ between cohorts (PFI, DFS, BCR), so the pooled estimate answers
     "does the module track progression in general", not one specific endpoint.

Outputs (in --output-dir): per_dataset.csv, meta.csv, analysis.md, heatmap.png, forest_top.png.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.module_scores import endpoint

Z95 = 1.959963984540054
MAX_COVARIATE_MISSING = 0.20      # drop a covariate if more than this share of eligible records lack it
MIN_EVENTS_PER_PARAMETER = 10     # below this the adjusted model is flagged as possibly overfit
MISSING_TOKENS = {'', 'nan', 'na', 'n/a', 'none', 'unknown', 'not available', 'x'}


# --------------------------------------------------------------------------- covariate parsing

def _missing(value) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value)) or str(value).strip().lower() in MISSING_TOKENS


def number(series: pd.Series) -> pd.Series:
    def parse(value):
        if _missing(value):
            return np.nan
        try:
            return float(value)
        except ValueError:
            raise ValueError(f'not a number: {value!r}') from None
    return series.map(parse).astype(float)


def grade_group(primary: int, secondary: int) -> float:
    """ISUP grade group from Gleason primary/secondary patterns."""
    total = primary + secondary
    if total <= 6:
        return 1.0
    if total == 7:
        return 2.0 if primary <= 3 else 3.0
    return 4.0 if total == 8 else 5.0


def grade_from_patterns(primary: pd.Series, secondary: pd.Series) -> pd.Series:
    p, s = number(primary), number(secondary)
    return pd.Series([np.nan if np.isnan(a) or np.isnan(b) else grade_group(int(a), int(b)) for a, b in zip(p, s)],
                     index=primary.index, dtype=float)


GLEASON_TEXT = re.compile(r'^\s*(\d+)\s*=\s*(\d)\s*\+\s*(\d)\s*$')


def grade_from_text(series: pd.Series) -> pd.Series:
    """Grade group from strings like '7=3+4' (Cambridge/Stockholm GEO characteristics)."""
    def parse(value):
        if _missing(value):
            return np.nan
        match = GLEASON_TEXT.match(str(value))
        if not match:
            raise ValueError(f'unrecognised Gleason value: {value!r}')
        total, primary, secondary = map(int, match.groups())
        if primary + secondary != total:
            raise ValueError(f'Gleason patterns do not add up: {value!r}')
        return grade_group(primary, secondary)
    return series.map(parse).astype(float)


T_STAGE = re.compile(r'^\s*P?T([0-4])', re.IGNORECASE)


def pt3_plus(series: pd.Series) -> pd.Series:
    """1 for pathological T3/T4, 0 for T0-T2 (T0: no residual tumour in the specimen), missing for Tx or blank."""
    def parse(value):
        if _missing(value) or re.match(r'^\s*P?TX', str(value), re.IGNORECASE):
            return np.nan
        match = T_STAGE.match(str(value))
        if not match:
            raise ValueError(f'unrecognised T stage: {value!r}')
        return 1.0 if int(match.group(1)) >= 3 else 0.0
    return series.map(parse).astype(float)


def gleason_8_plus(series: pd.Series) -> pd.Series:
    """1 for a total Gleason score of 8-10, 0 for 2-7. Used where only the total score is recorded."""
    scores = number(series)
    if ((scores < 2) | (scores > 10)).any():
        raise ValueError(f'Gleason score outside 2-10: {sorted(scores[(scores < 2) | (scores > 10)].unique())}')
    return (scores >= 8).astype(float).where(scores.notna())


def n1(series: pd.Series) -> pd.Series:
    def parse(value):
        if _missing(value) or str(value).strip().upper() == 'NX':
            return np.nan
        text = str(value).strip().upper()
        if text not in {'N0', 'N1'}:
            raise ValueError(f'unrecognised N stage: {value!r}')
        return float(text == 'N1')
    return series.map(parse).astype(float)


# Snapshot name prefix -> covariate -> (parser, source column(s)). Only verified columns are listed;
# a cohort that is not here is analysed unadjusted only.
COVARIATES = {
    'tcga-prad': {'age': (number, 'age'),
                  'grade_group': (grade_from_patterns, 'legacy_patient__GLEASON_PATTERN_PRIMARY',
                                  'legacy_patient__GLEASON_PATTERN_SECONDARY'),
                  'pt3_plus': (pt3_plus, 'path_t_group'),
                  'n1': (n1, 'patient__PATH_N_STAGE')},
    'taylor': {'grade_group': (grade_from_patterns, 'sample__GLEASON_SCORE_1', 'sample__GLEASON_SCORE_2'),
               'pt3_plus': (pt3_plus, 'patient__PATH_T_STAGE')},
    'cambridge': {'age': (number, 'geo__age_at_diag'),
                  'grade_group': (grade_from_text, 'geo__tumour_gleason'),
                  'pt3_plus': (pt3_plus, 'geo__pathology_stage')},
    'stockholm': {'grade_group': (grade_from_text, 'geo__tumour_gleason'),
                  'pt3_plus': (pt3_plus, 'geo__pathology_stage')},
    # Belfast biopsies: total Gleason only (no patterns, so no grade group) and clinical, not pathological, T stage.
    'belfast': {'age': (number, 'geo__patient_age_years'),
                'gleason_8_plus': (gleason_8_plus, 'geo__gleason_grade'),
                'ct3_plus': (pt3_plus, 'geo__t_stage')},
}


def covariate_spec(dataset: str) -> dict:
    return next((spec for prefix, spec in COVARIATES.items() if dataset.startswith(prefix)), {})


def clinical_covariates(clinical: pd.DataFrame, dataset: str) -> pd.DataFrame:
    frame = {}
    for name, (parser, *columns) in covariate_spec(dataset).items():
        absent = [c for c in columns if c not in clinical]
        if absent:
            raise ValueError(f'{dataset}: covariate {name} needs missing column(s) {absent}')
        try:
            frame[name] = parser(*(clinical[c] for c in columns))
        except ValueError as error:
            raise ValueError(f'{dataset}: covariate {name}: {error}') from None
    return pd.DataFrame(frame, index=clinical.index)


# --------------------------------------------------------------------------- statistics

def bh(p_values) -> np.ndarray:
    """Benjamini-Hochberg q-values; missing p-values stay missing."""
    p = np.asarray(p_values, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    n = int(ok.sum())
    if n:
        order = np.argsort(p[ok])
        ranked = p[ok][order] * n / np.arange(1, n + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        values = np.empty(n)
        values[order] = np.minimum(ranked, 1.0)
        q[ok] = values
    return q


def normal_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2))


def chi2_1_p(statistic: float) -> float:
    return math.erfc(math.sqrt(max(statistic, 0.0) / 2))


def random_effects(log_hr, se) -> dict:
    """DerSimonian-Laird pooled log hazard ratio across datasets."""
    y, v = np.asarray(log_hr, float), np.asarray(se, float) ** 2
    keep = np.isfinite(y) & np.isfinite(v) & (v > 0)
    y, v = y[keep], v[keep]
    k = len(y)
    if k < 2:
        return {'k': k}
    w = 1 / v
    fixed = (w * y).sum() / w.sum()
    q_stat = float((w * (y - fixed) ** 2).sum())
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (q_stat - (k - 1)) / c) if c > 0 else 0.0
    w_star = 1 / (v + tau2)
    pooled = float((w_star * y).sum() / w_star.sum())
    pooled_se = float(1 / math.sqrt(w_star.sum()))
    return {'k': k, 'hr': math.exp(pooled), 'ci_low': math.exp(pooled - Z95 * pooled_se),
            'ci_high': math.exp(pooled + Z95 * pooled_se), 'p': normal_p(pooled / pooled_se),
            'i2': max(0.0, (q_stat - (k - 1)) / q_stat) if q_stat > 0 else 0.0, 'tau2': tau2,
            'same_direction': int((np.sign(y) == np.sign(pooled)).sum())}


def cox(frame: pd.DataFrame):
    from lifelines import CoxPHFitter
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return CoxPHFitter().fit(frame, duration_col='months', event_col='event')


# --------------------------------------------------------------------------- per dataset

def adjusted_dataset(processed: Path, scores_path: Path, modules) -> tuple[list[dict], dict]:
    """Adjusted Cox rows for each module plus one record describing the dataset's model."""
    name = processed.name
    clinical = pd.read_csv(processed / 'clinical.csv', dtype={'sample_id': str}, low_memory=False)
    metadata = json.loads((processed / 'metadata.json').read_text())
    scores = pd.read_csv(scores_path, dtype={'sample_id': str}).set_index('sample_id')
    if scores.index.tolist() != clinical.sample_id.tolist():
        raise ValueError(f'{name}: {scores_path.name} samples do not match clinical.csv; rerun module-scores')
    endpoint_name, event, months, eligible = endpoint(clinical, metadata, name)
    covariates = clinical_covariates(clinical, name)
    covariates.index = clinical.sample_id
    eligible = pd.Series(eligible, index=clinical.sample_id)
    info = {'dataset': name, 'endpoint': endpoint_name, 'eligible': int(eligible.sum()),
            'covariates_dropped': ''}
    if covariates.empty:
        return [], {**info, 'covariates_used': '', 'note': 'no covariate mapping for this cohort; unadjusted only'}
    share_missing = covariates[eligible].isna().mean()
    dropped = [c for c in covariates if share_missing[c] > MAX_COVARIATE_MISSING or covariates.loc[eligible, c].nunique() < 2]
    covariates = covariates.drop(columns=dropped)
    info['covariates_dropped'] = ', '.join(dropped)
    base = pd.DataFrame({'months': months.to_numpy(float), 'event': event.to_numpy(float)}, index=clinical.sample_id)
    base = base.join(covariates)
    complete = eligible & base.notna().all(axis=1)
    base = base[complete]
    events = int(base.event.eq(1).sum())
    parameters = covariates.shape[1] + 1
    info.update({'covariates_used': ', '.join(covariates.columns), 'complete_cases': int(complete.sum()),
                 'events': events, 'events_per_parameter': round(events / parameters, 1),
                 'note': '' if events / parameters >= MIN_EVENTS_PER_PARAMETER
                 else f'fewer than {MIN_EVENTS_PER_PARAMETER} events per parameter; adjusted estimates may be unstable'})
    clinical_fit = cox(base)
    info['clinical_c_index'] = round(float(clinical_fit.concordance_index_), 4)
    rows = []
    for module in modules:
        if module not in scores:
            continue
        frame = base.assign(score=scores.loc[base.index, module].to_numpy(float))
        try:
            fit = cox(frame)
        except Exception as error:   # convergence failures are reported, not fatal
            rows.append({'dataset': name, 'module': module, 'adjusted_note': f'Cox fit failed: {error}'})
            continue
        row = fit.summary.loc['score']
        lr = 2 * (fit.log_likelihood_ - clinical_fit.log_likelihood_)
        rows.append({'dataset': name, 'module': module,
                     'adj_hazard_ratio_per_sd': float(row['exp(coef)']),
                     'adj_hr_ci_low': float(row['exp(coef) lower 95%']),
                     'adj_hr_ci_high': float(row['exp(coef) upper 95%']), 'adj_p_value': float(row['p']),
                     'added_value_lr_p': chi2_1_p(lr),
                     'delta_c_index': float(fit.concordance_index_) - info['clinical_c_index']})
    return rows, info


def per_dataset_table(summary: pd.DataFrame, adjusted: pd.DataFrame) -> pd.DataFrame:
    table = summary.copy()
    table['q_value'] = table.groupby('dataset').p_value.transform(lambda p: pd.Series(bh(p), index=p.index))
    if not adjusted.empty:
        table = table.merge(adjusted, on=['dataset', 'module'], how='left')
        table['adj_q_value'] = table.groupby('dataset').adj_p_value.transform(lambda p: pd.Series(bh(p), index=p.index))
    return table


def meta_table(table: pd.DataFrame, labels: dict) -> pd.DataFrame:
    rows = []
    for module, group in table.groupby('module'):
        record = {'module': module, 'label': labels.get(module, module)}
        for prefix, hr, low, high in (('', 'hazard_ratio_per_sd', 'hr_ci_low', 'hr_ci_high'),
                                      ('adj_', 'adj_hazard_ratio_per_sd', 'adj_hr_ci_low', 'adj_hr_ci_high')):
            if hr not in group:
                continue
            g = group.dropna(subset=[hr])
            result = random_effects(np.log(g[hr]), (np.log(g[high]) - np.log(g[low])) / (2 * Z95))
            record.update({f'{prefix}{key}': value for key, value in result.items()})
        rows.append(record)
    meta = pd.DataFrame(rows)
    for prefix in ('', 'adj_'):
        if f'{prefix}p' in meta:
            meta[f'{prefix}q'] = bh(meta[f'{prefix}p'])
    return meta.sort_values('p' if 'p' in meta else 'module').reset_index(drop=True)


# --------------------------------------------------------------------------- report and figures

SURFACE, INK, SECONDARY, MUTED, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
BLUE, MID, RED = '#2a78d6', '#f0efec', '#e34948'   # diverging pair: protective <-> higher risk


def fmt_hr(hr, low, high) -> str:
    return '' if pd.isna(hr) else f'{hr:.2f} ({low:.2f}–{high:.2f})'


def markdown(models: pd.DataFrame, table: pd.DataFrame, meta: pd.DataFrame, top: int, excluded) -> str:
    lines = ['# Module analysis', '',
             f'Datasets: {", ".join(sorted(table.dataset.unique()))}'
             + (f' (excluded from meta-analysis: {", ".join(excluded)})' if excluded else ''), '',
             '## Clinical models', '',
             '| Dataset | Endpoint | Covariates used | Dropped | Complete cases (events) | Events/parameter | Clinical C-index | Note |',
             '|---|---|---|---|---:|---:|---:|---|']
    for r in models.itertuples():
        lines.append(f'| {r.dataset} | {str(r.endpoint).upper()} | {r.covariates_used or "—"} | '
                     f'{r.covariates_dropped or "—"} | {getattr(r, "complete_cases", "")} ({getattr(r, "events", "")}) | '
                     f'{getattr(r, "events_per_parameter", "")} | {getattr(r, "clinical_c_index", "")} | {r.note} |')
    lines += ['', '## Significant modules per dataset (BH q < 0.05)', '', '| Dataset | Scored | Unadjusted | Adjusted |',
              '|---|---:|---:|---:|']
    for dataset, g in table.dropna(subset=['p_value']).groupby('dataset'):
        adjusted = int((g.adj_q_value < 0.05).sum()) if 'adj_q_value' in g else 0
        lines.append(f'| {dataset} | {len(g)} | {int((g.q_value < 0.05).sum())} | {adjusted} |')
    lines += ['', f'## Top {top} modules, random-effects meta-analysis', '',
              '| Module | k | Pooled HR per SD (95% CI) | q | I² | Same direction | Adjusted pooled HR (95% CI) | Adjusted q |',
              '|---|---:|---|---:|---:|---:|---|---:|']
    for r in meta.head(top).itertuples():
        adjusted = fmt_hr(getattr(r, 'adj_hr', np.nan), getattr(r, 'adj_ci_low', np.nan), getattr(r, 'adj_ci_high', np.nan))
        adj_q = getattr(r, 'adj_q', np.nan)
        lines.append(f'| {r.label} | {int(r.k)} | {fmt_hr(r.hr, r.ci_low, r.ci_high)} | {r.q:.2g} | {r.i2:.0%} | '
                     f'{int(r.same_direction)}/{int(r.k)} | {adjusted} | {"" if pd.isna(adj_q) else f"{adj_q:.2g}"} |')
    lines += ['', 'HR > 1: higher score, higher risk. Endpoints differ by cohort and are not relabeled.',
              'Adjusted models use complete cases and the covariates listed above; C-indices are in-sample.']
    return '\n'.join(lines) + '\n'


def _style(ax):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=SECONDARY, labelsize=8, length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def heatmap(table: pd.DataFrame, meta: pd.DataFrame, top: int, path: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    modules = meta.dropna(subset=['p']).head(top)
    datasets = sorted(table.dataset.unique())
    grid = table.pivot(index='module', columns='dataset', values='hazard_ratio_per_sd').reindex(modules.module)[datasets]
    q = table.pivot(index='module', columns='dataset', values='q_value').reindex(modules.module)[datasets]
    values = np.log(grid.to_numpy(float))
    limit = float(np.nanmax(np.abs(values))) if np.isfinite(values).any() else 1.0
    cmap = LinearSegmentedColormap.from_list('hr', [BLUE, MID, RED]).with_extremes(bad=SURFACE)
    fig, ax = plt.subplots(figsize=(1.5 * len(datasets) + 3.2, 0.34 * len(modules) + 1.8), facecolor=SURFACE)
    _style(ax)
    mesh = ax.pcolormesh(np.ma.masked_invalid(values), cmap=cmap, norm=TwoSlopeNorm(0, -limit, limit),
                         edgecolors=SURFACE, linewidth=2)
    for i in range(len(modules)):
        for j in range(len(datasets)):
            hr = grid.iat[i, j]
            if pd.notna(hr):
                star = '*' if q.iat[i, j] < 0.05 else ''
                ax.text(j + 0.5, i + 0.5, f'{hr:.2f}{star}', ha='center', va='center', fontsize=7.5, color=INK)
    ax.set_xticks(np.arange(len(datasets)) + 0.5, datasets)
    ax.set_yticks(np.arange(len(modules)) + 0.5, modules.label)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    bar = fig.colorbar(mesh, ax=ax, fraction=0.04, pad=0.02)
    ticks = [t for t in (0.5, 0.7, 1, 1.4, 2) if abs(math.log(t)) <= limit]
    bar.set_ticks([math.log(t) for t in ticks], labels=[f'{t:g}' for t in ticks])
    bar.outline.set_visible(False)
    bar.ax.tick_params(colors=SECONDARY, labelsize=8, length=0)
    bar.set_label('HR per SD (blue = lower risk, red = higher risk)', color=SECONDARY, fontsize=8)
    fig.suptitle(f'Top {len(modules)} modules by meta-analysis p; cell = unadjusted HR, * = BH q < 0.05',
                 color=INK, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def forest(table: pd.DataFrame, meta: pd.DataFrame, top: int, path: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    modules = meta.dropna(subset=['p']).head(top)
    datasets = sorted(table.dataset.unique(), reverse=True)
    rows = ['pooled'] + datasets
    cols = min(4, len(modules))
    nrows = math.ceil(len(modules) / cols)
    fig, axes = plt.subplots(nrows, cols, figsize=(3.0 * cols + 1.0, (0.32 * len(rows) + 1.0) * nrows),
                             sharex=True, sharey=True, facecolor=SURFACE, squeeze=False)
    low = min(table.hr_ci_low.min(), modules.ci_low.min(), 0.5)
    high = max(table.hr_ci_high.max(), modules.ci_high.max(), 2.0)
    for ax, record in zip(axes.flat, modules.itertuples()):
        _style(ax)
        ax.set_xscale('log')
        ax.set_xlim(low / 1.1, high * 1.1)
        ax.axvline(1.0, color=MUTED, linewidth=1, zorder=1)
        ax.grid(axis='x', color=GRID, linewidth=0.6, zorder=0)
        for y, dataset in enumerate(rows):
            if dataset == 'pooled':
                ax.plot([record.ci_low, record.ci_high], [y, y], color=INK, linewidth=2, zorder=2)
                ax.plot(record.hr, y, 'D', markersize=8, color=INK, zorder=3)
                continue
            hit = table[(table.module == record.module) & (table.dataset == dataset)].dropna(subset=['hazard_ratio_per_sd'])
            if hit.empty:
                continue
            r = hit.iloc[0]
            ax.plot([r.hr_ci_low, r.hr_ci_high], [y, y], color=BLUE, linewidth=2, solid_capstyle='round', zorder=2)
            ax.plot(r.hazard_ratio_per_sd, y, 'o', markersize=8, markeredgewidth=2, markeredgecolor=BLUE,
                    markerfacecolor=BLUE if r.q_value < 0.05 else SURFACE, zorder=3)
        ax.set_title(record.label, color=INK, fontsize=9)
        ax.set_yticks(range(len(rows)), rows)
        ticks = [t for t in (0.25, 0.5, 1, 2, 4) if low / 1.1 <= t <= high * 1.1]
        ax.set_xticks(ticks, [f'{t:g}' for t in ticks])
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for ax in list(axes.flat)[len(modules):]:
        ax.set_visible(False)
    fig.supxlabel('HR per SD of module score (log scale). Circles: datasets, filled = BH q < 0.05. '
                  'Diamond: random-effects pooled.', color=SECONDARY, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


# --------------------------------------------------------------------------- command line

def load_labels(path: Path | None) -> dict:
    if path is None or not path.is_file():
        return {}
    info = pd.read_csv(path, dtype=str).fillna('')
    return {r.module: f'{r.module} · {r.annotation}' if r.annotation else r.module for r in info.itertuples()}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--scores-dir', type=Path, help='module-scores output; default <data-dir>/analysis/module_scores')
    parser.add_argument('--module-info', type=Path, help='module,category,annotation CSV for readable labels')
    parser.add_argument('--exclude', action='append', default=[], help='Leave this snapshot out of the meta-analysis')
    parser.add_argument('--top', type=int, default=20, help='Modules shown in the report and heatmap')
    parser.add_argument('--output-dir', type=Path, help='Default: <scores-dir>/analysis')
    args = parser.parse_args(argv)
    scores_dir = args.scores_dir or args.data_dir / 'analysis' / 'module_scores'
    summary_path = scores_dir / 'summary.csv'
    if not summary_path.is_file():
        parser.error(f'{summary_path} not found; run module-scores first')
    summary = pd.read_csv(summary_path)
    output = args.output_dir or scores_dir / 'analysis'
    output.mkdir(parents=True, exist_ok=True)
    models, adjusted = [], []
    for dataset, group in summary.groupby('dataset'):
        rows, info = adjusted_dataset(args.data_dir / 'processed' / dataset, scores_dir / f'{dataset}_scores.csv',
                                      group.dropna(subset=['p_value']).module.tolist())
        adjusted += rows
        models.append(info)
    models = pd.DataFrame(models)
    table = per_dataset_table(summary, pd.DataFrame(adjusted))
    labels = load_labels(args.module_info)
    meta = meta_table(table[~table.dataset.isin(args.exclude)].dropna(subset=['p_value']), labels)
    table.to_csv(output / 'per_dataset.csv', index=False, float_format='%.6g')
    meta.to_csv(output / 'meta.csv', index=False, float_format='%.6g')
    models.to_csv(output / 'clinical_models.csv', index=False)
    shown = table.dropna(subset=['p_value'])
    heatmap(shown, meta, args.top, output / 'heatmap.png')
    forest(shown[~shown.dataset.isin(args.exclude)], meta, min(args.top, 8), output / 'forest_top.png')
    report = markdown(models, table, meta, args.top, args.exclude)
    (output / 'analysis.md').write_text(report, encoding='utf-8')
    print(report, end='')
    print(f'\nWrote {output}/per_dataset.csv, meta.csv, clinical_models.csv, analysis.md, heatmap.png, forest_top.png')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
