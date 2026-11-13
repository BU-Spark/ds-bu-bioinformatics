"""PCA gene-module scores for each processed dataset, tested against that dataset's own endpoint.

For every dataset and module:
  1. Match module genes to genes.csv gene_symbol. When several feature rows (probes) share a symbol,
     use the row with the highest mean expression, so a gene with many probes does not dominate.
     Rows with no variance are dropped. The output files are not changed; this happens only here.
  2. Expression that is not on a log scale (max > 50, currently Taylor) is log2(x + 1) transformed
     here, as a placeholder until the team decides Taylor's transform.
  3. Standardize each gene across samples, then take the first principal component (PC1). Its sign
     is set so that a higher score means higher average expression of the module genes. Scores are
     scaled to mean 0, SD 1 within the dataset.
  4. Fit a Cox model of the dataset's endpoint on the score (hazard ratio per SD), using eligible
     records only, and report Harrell's concordance index.

Endpoints are not relabeled: TCGA uses PFI, Taylor DFS, Cambridge and Stockholm BCR, and Sboner
prostate-cancer-specific mortality. Scores are therefore compared within each dataset, not pooled.
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

LOG_SCALE_MAX = 50.0
MIN_GENES = 5
# Sboner (William's branch) stores one follow-up time and eligibility flag for all of its endpoints
# instead of <endpoint>_months/_eligible. Used only until that output adopts the shared pattern.
SHARED_FOLLOWUP = ('followup_months', 'survival_eligible')
FALLBACK_ENDPOINT = {'sboner': 'pcsm'}


def load_modules(path: Path) -> dict[str, list[str]]:
    frame = pd.read_csv(path, dtype=str)
    if not {'gene', 'module'}.issubset(frame.columns):
        raise ValueError('Module file needs gene and module columns')
    frame['gene'] = frame.gene.str.strip()
    frame = frame.dropna().drop_duplicates()
    return {module: group.gene.tolist() for module, group in frame.groupby('module', sort=True)}


def load_dataset(processed: Path):
    expression = pd.read_csv(processed / 'expression.csv.gz', index_col=0, dtype={'feature_id': str})
    expression.index = expression.index.astype(str)
    genes = pd.read_csv(processed / 'genes.csv', dtype={'feature_id': str, 'gene_id': str, 'gene_symbol': str})
    clinical = pd.read_csv(processed / 'clinical.csv', dtype={'sample_id': str}, low_memory=False)
    metadata = json.loads((processed / 'metadata.json').read_text())
    if genes.feature_id.tolist() != expression.index.tolist() or clinical.sample_id.tolist() != expression.columns.tolist():
        raise ValueError(f'{processed.name}: expression, genes.csv and clinical.csv are not aligned; run check-outputs')
    transform = 'none'
    if float(expression.to_numpy().max()) > LOG_SCALE_MAX:
        expression = np.log2(expression.clip(lower=0) + 1)
        transform = 'log2(x + 1) applied for scoring (deposited values are not log-scaled)'
    return expression, genes, clinical, metadata, transform


def module_matrix(expression: pd.DataFrame, genes: pd.DataFrame, symbols) -> pd.DataFrame:
    """Genes x samples for one module: one feature per symbol (highest mean), non-constant rows only."""
    wanted = genes.loc[genes.gene_symbol.isin(set(symbols)), ['feature_id', 'gene_symbol']]
    if wanted.empty:
        return pd.DataFrame(columns=expression.columns)
    rows = expression.loc[wanted.feature_id]
    order = pd.DataFrame({'symbol': wanted.gene_symbol.to_numpy(), 'feature_id': wanted.feature_id.to_numpy(),
                          'mean': rows.mean(axis=1).to_numpy()})
    best = order.sort_values(['symbol', 'mean', 'feature_id'], ascending=[True, False, True]).drop_duplicates('symbol')
    matrix = expression.loc[best.feature_id]
    matrix.index = pd.Index(best.symbol.tolist(), name='gene_symbol')
    return matrix[matrix.std(axis=1, ddof=1) > 0]


def pca_score(matrix: pd.DataFrame):
    """PC1 of the standardized genes x samples matrix. Returns (score per sample, variance share, loadings)."""
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=1), axis=0)
    _, singular, vt = np.linalg.svd(z.to_numpy().T, full_matrices=False)
    loadings = pd.Series(vt[0], index=matrix.index)
    score = pd.Series(z.to_numpy().T @ vt[0], index=matrix.columns)
    if np.corrcoef(score, z.mean(axis=0))[0, 1] < 0:   # higher score = higher module expression
        score, loadings = -score, -loadings
    variance = float(singular[0] ** 2 / (singular ** 2).sum())
    return (score - score.mean()) / score.std(ddof=1), variance, loadings


def endpoint(clinical: pd.DataFrame, metadata: dict, dataset: str):
    """Return (name, event, months, eligible) for the dataset's primary endpoint."""
    prefixes = [c[:-6] for c in clinical.columns if c.endswith('_event') and '__' not in c
                and f'{c[:-6]}_months' in clinical and f'{c[:-6]}_eligible' in clinical]
    name = metadata.get('primary_endpoint') if metadata.get('primary_endpoint') in prefixes else None
    if name is None and len(prefixes) == 1:
        name = prefixes[0]
    if name is not None:
        months, eligible = clinical[f'{name}_months'], clinical[f'{name}_eligible']
    else:
        name = next((v for k, v in FALLBACK_ENDPOINT.items() if dataset.startswith(k)), None)
        if name is None or f'{name}_event' not in clinical or not set(SHARED_FOLLOWUP).issubset(clinical.columns):
            raise ValueError(f'{dataset}: cannot identify a primary endpoint (found {prefixes})')
        months, eligible = clinical[SHARED_FOLLOWUP[0]], clinical[SHARED_FOLLOWUP[1]]
    eligible = eligible.map(lambda v: v in {True, 'True', 'true', 1})
    return (name, pd.to_numeric(clinical[f'{name}_event'], errors='coerce'),
            pd.to_numeric(months, errors='coerce'), eligible.to_numpy(dtype=bool))


def survival(score: pd.Series, months: pd.Series, event: pd.Series) -> dict:
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    frame = pd.DataFrame({'score': score.to_numpy(), 'months': months.to_numpy(dtype=float),
                          'event': event.to_numpy(dtype=float)})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fit = CoxPHFitter().fit(frame, duration_col='months', event_col='event')
    row = fit.summary.loc['score']
    return {'hazard_ratio_per_sd': float(row['exp(coef)']), 'hr_ci_low': float(row['exp(coef) lower 95%']),
            'hr_ci_high': float(row['exp(coef) upper 95%']), 'p_value': float(row['p']),
            # A higher score predicting earlier events gives a concordance above 0.5.
            'concordance': float(concordance_index(frame.months, -frame.score, frame.event))}


def score_dataset(processed: Path, modules: dict[str, list[str]]):
    expression, genes, clinical, metadata, transform = load_dataset(processed)
    name, event, months, eligible = endpoint(clinical, metadata, processed.name)
    scores, rows = {}, []
    for module, symbols in modules.items():
        matrix = module_matrix(expression, genes, symbols)
        row = {'dataset': processed.name, 'module': module, 'module_genes': len(symbols), 'genes_used': len(matrix),
               'endpoint': name, 'eligible': int(eligible.sum()), 'events': int(event[eligible].eq(1).sum()),
               'expression_transform': transform}
        if len(matrix) < MIN_GENES:
            rows.append({**row, 'note': f'fewer than {MIN_GENES} module genes measured; not scored'})
            continue
        score, variance, _ = pca_score(matrix)
        scores[module] = score
        row['pc1_variance_explained'] = round(variance, 4)
        row.update(survival(score[eligible], months[eligible], event[eligible]))
        rows.append(row)
    return pd.DataFrame(scores).rename_axis('sample_id'), pd.DataFrame(rows)


def markdown(summary: pd.DataFrame) -> str:
    lines = ['| Dataset | Endpoint (eligible, events) | Module | Genes used | PC1 variance | HR per SD (95% CI) | p | C-index |',
             '|---|---|---|---:|---:|---|---:|---:|']
    for r in summary.itertuples():
        if pd.isna(getattr(r, 'hazard_ratio_per_sd', np.nan)):
            lines.append(f'| {r.dataset} | {r.endpoint.upper()} ({r.eligible}, {r.events}) | {r.module} | '
                         f'{r.genes_used}/{r.module_genes} | | not scored | | |')
            continue
        lines.append(f'| {r.dataset} | {r.endpoint.upper()} ({r.eligible}, {r.events}) | {r.module} | '
                     f'{r.genes_used}/{r.module_genes} | {r.pc1_variance_explained:.0%} | '
                     f'{r.hazard_ratio_per_sd:.2f} ({r.hr_ci_low:.2f}–{r.hr_ci_high:.2f}) | {r.p_value:.2g} | '
                     f'{r.concordance:.2f} |')
    return '\n'.join(lines) + '\n'


def forest_plot(summary: pd.DataFrame, path: Path) -> None:
    """Hazard ratio per SD with 95% CI: one panel per module, one row per dataset, log scale."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    surface, ink, secondary, muted, grid, blue = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#2a78d6'
    data = summary.dropna(subset=['hazard_ratio_per_sd'])
    modules = sorted(data.module.unique())
    datasets = sorted(data.dataset.unique(), reverse=True)
    low = min(data.hr_ci_low.min(), 0.5)
    high = max(data.hr_ci_high.max(), 2.0)
    fig, axes = plt.subplots(1, len(modules), figsize=(2.6 * len(modules) + 1.4, 0.45 * len(datasets) + 1.6),
                             sharey=True, facecolor=surface)
    axes = np.atleast_1d(axes)
    for ax, module in zip(axes, modules):
        ax.set_facecolor(surface)
        ax.set_xscale('log')
        ax.set_xlim(low / 1.15, high * 1.15)
        ax.axvline(1.0, color=muted, linewidth=1, zorder=1)
        ax.grid(axis='x', color=grid, linewidth=0.6, zorder=0)
        for y, dataset in enumerate(datasets):
            hit = data[(data.module == module) & (data.dataset == dataset)]
            if hit.empty:
                continue
            r = hit.iloc[0]
            ax.plot([r.hr_ci_low, r.hr_ci_high], [y, y], color=blue, linewidth=2, solid_capstyle='round', zorder=2)
            significant = r.p_value < 0.05
            ax.plot(r.hazard_ratio_per_sd, y, 'o', markersize=8, zorder=3, markeredgewidth=2,
                    markeredgecolor=blue, markerfacecolor=blue if significant else surface)
        ax.set_title(module.replace('_', ' '), color=ink, fontsize=10)
        ax.tick_params(colors=secondary, labelsize=8, length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_yticks(range(len(datasets)), datasets)
        ticks = [t for t in (0.25, 0.5, 1, 2, 4) if low / 1.15 <= t <= high * 1.15]
        ax.set_xticks(ticks, [f'{t:g}' for t in ticks])
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    fig.supxlabel('Hazard ratio per SD of module score (log scale; filled = p < 0.05)', color=secondary, fontsize=9)
    fig.suptitle('PCA module score vs. each dataset\'s endpoint (placeholder Hallmark modules)', color=ink, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=surface)
    plt.close(fig)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--modules', type=Path, default=Path('modules/hallmark_placeholder.csv'))
    parser.add_argument('--snapshot', action='append', help='Score only this snapshot (repeatable); default: all')
    parser.add_argument('--output-dir', type=Path, help='Default: <data-dir>/analysis/module_scores')
    args = parser.parse_args(argv)
    modules = load_modules(args.modules)
    root = args.data_dir / 'processed'
    names = args.snapshot or sorted(p.name for p in root.iterdir() if p.is_dir() and not p.name.startswith('.'))
    output = args.output_dir or args.data_dir / 'analysis' / 'module_scores'
    output.mkdir(parents=True, exist_ok=True)
    summaries = []
    for name in names:
        scores, summary = score_dataset(root / name, modules)
        scores.to_csv(output / f'{name}_scores.csv', float_format='%.6g')
        summaries.append(summary)
    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(output / 'summary.csv', index=False, float_format='%.6g')
    forest_plot(summary, output / 'hazard_ratios.png')
    (output / 'summary.md').write_text(markdown(summary))
    print(markdown(summary), end='')
    print(f'\nWrote {output}/summary.csv, summary.md, hazard_ratios.png and <snapshot>_scores.csv')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
