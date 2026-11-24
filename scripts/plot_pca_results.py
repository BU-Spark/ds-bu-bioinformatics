"""Figures for the PCA module scores written by pca_all_cohorts.

Reads <pca-dir>/coverage.csv, cohorts.csv and <cohort>_scores.csv (outputs only; nothing is recomputed
except the PC1/PC2 view in figure 1) and writes four PNGs to <pca-dir>/figures/:

  1_pca_explainer.png        PC1 vs PC2 for one module in one cohort (what a module score is)
  2_module_coverage.png      share of each module's genes measured, cohorts x modules
  3_pc1_variance.png         variance captured by PC1 for every scored module, per cohort
  4_scores_heatmap.png       module scores per patient for selected modules in one cohort

    python scripts/plot_pca_results.py --pca-dir data/analysis/pca_all --data-dir data/full-validation
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from ds_bu_bioinformatics.module_scores import load_dataset, load_modules, module_matrix  # noqa: E402

SURF, INK, SEC, MUTED, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
BLUE, MID, RED = '#2a78d6', '#f0efec', '#e34948'
FOCUS = ['L12', 'L11', 'L147', 'L127', 'L58', 'L59', 'L72', 'L16', 'L60']


def style(ax, grid_axis='y'):
    ax.set_facecolor(SURF)
    ax.tick_params(colors=SEC, labelsize=8, length=0)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.6, zorder=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def colorbar(fig, mesh, ax, text):
    bar = fig.colorbar(mesh, ax=ax, fraction=0.02, pad=0.01)
    bar.set_label(text, color=SEC, fontsize=8)
    bar.outline.set_visible(False)
    bar.ax.tick_params(colors=SEC, labelsize=8, length=0)


def short_name(dataset: str) -> str:
    return dataset.replace('-eopc-public-recovery', '').rsplit('-v', 1)[0]


def explainer(processed: Path, modules_path: Path, module: str, label, path: Path) -> bool:
    if not processed.is_dir():
        return False
    expr, genes, clinical, _, _ = load_dataset(processed)
    matrix = module_matrix(expr, genes, load_modules(modules_path).get(module, []))
    if len(matrix) < 3:
        return False
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=1), axis=0).to_numpy().T
    _, singular, vt = np.linalg.svd(z, full_matrices=False)
    pcs = z @ vt[:2].T
    if np.corrcoef(pcs[:, 0], z.mean(axis=1))[0, 1] < 0:   # same orientation rule as module-scores
        pcs[:, 0] *= -1
    share = singular ** 2 / (singular ** 2).sum()
    fig, ax = plt.subplots(figsize=(7.2, 5.2), facecolor=SURF)
    style(ax, None)
    gleason = next((pd.to_numeric(clinical[c], errors='coerce').to_numpy() for c in
                    ('gleason_score', 'legacy_patient__GLEASON_SCORE', 'sample__GLEASON_SCORE') if c in clinical), None)
    if gleason is not None and np.isfinite(gleason).any():
        for mask, color, name in ((gleason <= 6, '#9fc3ee', 'Gleason ≤6'), (gleason == 7, '#5b9ae3', 'Gleason 7'),
                                  (gleason >= 8, '#1b5fb0', 'Gleason 8–10')):
            if mask.any():
                ax.scatter(pcs[mask, 0], pcs[mask, 1], s=18, color=color, edgecolor=SURF, linewidth=0.6,
                           label=f'{name} (n={int(mask.sum())})', zorder=3)
        ax.legend(frameon=False, fontsize=8, loc='upper left')
    else:
        ax.scatter(pcs[:, 0], pcs[:, 1], s=18, color=BLUE, edgecolor=SURF, linewidth=0.6, zorder=3)
    ax.axhline(0, color=GRID, linewidth=0.8, zorder=1)
    ax.axvline(0, color=GRID, linewidth=0.8, zorder=1)
    ax.set_xlabel(f'PC1 = module score ({share[0]:.0%} of variance) → higher module activity', color=SEC, fontsize=9)
    ax.set_ylabel(f'PC2 ({share[1]:.0%} of variance)', color=SEC, fontsize=9)
    ax.set_title(f'What PCA does: {matrix.shape[0]} genes of {label(module)} → one score per patient\n'
                 f'{short_name(processed.name)}, {matrix.shape[1]} tumors; each dot is a patient', color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)
    return True


def coverage_heatmap(cov: pd.DataFrame, cohorts: pd.DataFrame, path: Path) -> None:
    cov = cov.assign(share=cov.genes_used / cov.module_genes)
    modules = sorted(cov.module.unique(), key=lambda m: int(m[1:]) if m[1:].isdigit() else 0)
    order = cohorts.sort_values('modules_scored', ascending=False).dataset
    share = cov.pivot(index='dataset', columns='module', values='share').reindex(index=order, columns=modules)
    scored = cov.pivot(index='dataset', columns='module', values='scored').reindex(index=order, columns=modules)
    counts = cohorts.set_index('dataset').modules_scored
    fig, ax = plt.subplots(figsize=(13, 0.32 * len(order) + 1.4), facecolor=SURF)
    style(ax, None)
    cmap = LinearSegmentedColormap.from_list('cov', ['#f0efec', '#9fc3ee', BLUE, '#0d3b73']).with_extremes(bad='#ffffff')
    values = np.ma.masked_where(~scored.fillna(False).to_numpy(dtype=bool), share.to_numpy(dtype=float))
    mesh = ax.pcolormesh(values, cmap=cmap, vmin=0, vmax=1, edgecolors=SURF, linewidth=0.3)
    ax.set_yticks(np.arange(len(order)) + 0.5, [f'{short_name(d)} ({int(counts[d])})' for d in order])
    ticks = [i for i, m in enumerate(modules) if m[1:].isdigit() and int(m[1:]) % 10 == 0]
    ax.set_xticks(np.array(ticks) + 0.5, [modules[i] for i in ticks])
    ax.invert_yaxis()
    colorbar(fig, mesh, ax, 'Share of module genes measured')
    ax.set_title('Module coverage across cohorts (white = not scored, too few genes; modules scored in brackets)',
                 color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)


def variance_plot(cov: pd.DataFrame, path: Path) -> None:
    scored = cov[cov.scored.astype(bool)]
    order = scored.groupby('dataset').pc1_variance_explained.median().sort_values().index
    data = [scored.loc[scored.dataset == d, 'pc1_variance_explained'].to_numpy() * 100 for d in order]
    fig, ax = plt.subplots(figsize=(8.5, 0.37 * len(order) + 1.6), facecolor=SURF)
    style(ax, 'x')
    ax.boxplot(data, orientation='horizontal', widths=0.55, patch_artist=True, showfliers=False,
               medianprops=dict(color=INK, linewidth=1.5), whiskerprops=dict(color=MUTED),
               capprops=dict(color=MUTED), boxprops=dict(facecolor='#b9d4f3', edgecolor=BLUE))
    for i, values in enumerate(data, start=1):
        jitter = np.random.default_rng(i).uniform(-0.18, 0.18, len(values))
        ax.scatter(values, i + jitter, s=6, color=BLUE, alpha=0.45, zorder=3, linewidth=0)
    ax.set_yticks(range(1, len(order) + 1), [short_name(d) for d in order])
    ax.set_xlabel('Variance of module genes captured by PC1 (%), each dot is a module', color=SEC, fontsize=9)
    ax.set_title('How well one PCA score summarises each module, per cohort', color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)


def scores_heatmap(scores_path: Path, label, path: Path, focus=FOCUS) -> bool:
    if not scores_path.is_file():
        return False
    scores = pd.read_csv(scores_path, index_col=0)
    focus = [m for m in focus if m in scores]
    if not focus:
        return False
    scores = scores[focus].sort_values(focus[0])
    fig, ax = plt.subplots(figsize=(12, 0.4 * len(focus) + 1.2), facecolor=SURF)
    style(ax, None)
    cmap = LinearSegmentedColormap.from_list('sc', [BLUE, MID, RED])
    mesh = ax.pcolormesh(scores.T.to_numpy(), cmap=cmap, norm=TwoSlopeNorm(0, -2.5, 2.5))
    ax.set_yticks(np.arange(len(focus)) + 0.5, [label(m) for m in focus])
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_xlabel(f'{len(scores)} patients in {short_name(scores_path.name.replace("_scores.csv", ""))}, '
                  f'ordered by {label(focus[0])} →', color=SEC, fontsize=9)
    colorbar(fig, mesh, ax, 'Module score (SD)')
    ax.set_title('PCA module scores per patient (red = module more active, blue = less active)', color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor=SURF)
    plt.close(fig)
    return True


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--pca-dir', type=Path, required=True, help='pca_all_cohorts output folder')
    parser.add_argument('--data-dir', type=Path, default=Path('data/full-validation'), help='Root with processed/<snapshot>/')
    parser.add_argument('--modules', type=Path, default=Path('modules/camplab_modules.csv'))
    parser.add_argument('--module-info', type=Path, default=Path('modules/camplab_modules_info.csv'))
    parser.add_argument('--cohort', default='tcga-prad-v2', help='Cohort for figures 1 and 4')
    parser.add_argument('--module', default='L12', help='Module for figure 1')
    parser.add_argument('--output-dir', type=Path, help='Default: <pca-dir>/figures')
    args = parser.parse_args(argv)
    cov_path, cohorts_path = args.pca_dir / 'coverage.csv', args.pca_dir / 'cohorts.csv'
    if not cov_path.is_file() or not cohorts_path.is_file():
        parser.error(f'{args.pca_dir} has no coverage.csv/cohorts.csv; run pca_all_cohorts first')
    names = {}
    if args.module_info.is_file():
        info = pd.read_csv(args.module_info, dtype=str).fillna('')
        names = dict(zip(info.module, info.annotation.str.replace('Metallotheoin', 'Metallothionein')))
    label = lambda m: f'{m} · {names[m]}' if names.get(m) else m
    out = args.output_dir or args.pca_dir / 'figures'
    out.mkdir(parents=True, exist_ok=True)
    cov, cohorts = pd.read_csv(cov_path), pd.read_csv(cohorts_path)
    written = []
    if explainer(args.data_dir / 'processed' / args.cohort, args.modules, args.module, label, out / '1_pca_explainer.png'):
        written.append('1_pca_explainer.png')
    else:
        print(f'Skipped figure 1: {args.cohort} or module {args.module} not available')
    coverage_heatmap(cov, cohorts, out / '2_module_coverage.png')
    variance_plot(cov, out / '3_pc1_variance.png')
    written += ['2_module_coverage.png', '3_pc1_variance.png']
    if scores_heatmap(args.pca_dir / f'{args.cohort}_scores.csv', label, out / '4_scores_heatmap.png'):
        written.append('4_scores_heatmap.png')
    else:
        print(f'Skipped figure 4: no scores for {args.cohort}')
    print(f'Wrote {", ".join(written)} to {out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
