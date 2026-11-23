"""Within-cohort PCA versus GSVA on complete, structurally checked exports."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .module_scores import MIN_GENES, load_dataset, load_modules, module_matrix, pca_score
from .output_check import check_dataset
from .provenance import code_hashes
from .sources import sha256, write_json


def score_expression(expression, genes, modules, *, threads=1):
    """One outcome-independent representative per symbol; common variable-gene universe.

    PCA standardizes each module gene; GSVA uses continuous log-expression
    against the full measured variable-gene background, never module-only data.
    """
    import gseapy

    if expression.shape[1] < 3 or not np.isfinite(expression.to_numpy()).all():
        raise ValueError('Need at least three samples and complete finite expression')
    universe = module_matrix(expression, genes, genes.gene_symbol.dropna().unique())
    if universe.empty or not universe.index.is_unique:
        raise ValueError('No unique variable gene-symbol background')
    measured = set(universe.index)
    sets, pca, stats, loadings = {}, {}, [], []
    for module, symbols in modules.items():
        matched = sorted(set(symbols) & measured)
        row = {'module': module, 'module_genes': len(set(symbols)), 'genes_used': len(matched),
               'coverage_fraction': len(matched) / len(set(symbols)), 'status': 'scored'}
        if len(matched) < MIN_GENES or len(matched) == len(universe):
            row['status'] = 'not scored: fewer than 5 variable genes or no background genes'
            stats.append(row)
            continue
        matrix = universe.loc[matched]
        score, variance, weights = pca_score(matrix)
        z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=1), axis=0)
        average = z.mean(axis=0)
        orientation = float(np.corrcoef(score, average)[0, 1]) if average.std() > 1e-12 else np.nan
        row.update(pc1_variance_explained=variance, pc1_correlation_with_mean_gene_z=orientation,
                   pc1_orientation='higher mean gene expression' if np.isfinite(orientation) else 'undefined: mean gene-z is constant',
                   positive_loading_fraction=float(weights.gt(0).mean()))
        pca[module], sets[module] = score, matched
        loadings.extend({'module': module, 'gene_symbol': gene, 'pc1_loading': float(value)}
                        for gene, value in weights.items())
        stats.append(row)
    pca = pd.DataFrame(pca, index=expression.columns).rename_axis('sample_id')
    if not sets:
        return pca, pca.copy(), pd.DataFrame(stats), pd.DataFrame(loadings), universe
    result = gseapy.gsva(data=universe, gene_sets=sets, kcdf='Gaussian', weight=1.0,
                        mx_diff=True, abs_rnk=False, min_size=MIN_GENES,
                        max_size=max(map(len, sets.values())), threads=threads, seed=123,
                        outdir=None, verbose=False)
    gsva = result.res2d.pivot(index='Name', columns='Term', values='ES')
    if set(gsva.index) != set(pca.index) or set(gsva.columns) != set(pca.columns):
        raise ValueError('GSVA silently dropped/added samples or accepted gene sets')
    gsva = gsva.reindex(index=pca.index, columns=pca.columns).astype(float).rename_axis('sample_id')
    gsva.columns.name = None
    if not np.isfinite(gsva.to_numpy()).all() or not np.isfinite(pca.to_numpy()).all():
        raise ValueError('Nonfinite PCA/GSVA scores')
    for row in stats:
        module = row['module']
        if module not in pca:
            continue
        if pca[module].nunique() < 2 or gsva[module].nunique() < 2:
            row['comparison_note'] = 'Undefined correlation: constant scores'
            continue
        row['pearson_r'] = float(pca[module].corr(gsva[module]))
        row['spearman_rho'] = float(spearmanr(pca[module], gsva[module]).statistic)
    return pca, gsva, pd.DataFrame(stats), pd.DataFrame(loadings), universe


def figures(summary, output, scores):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    modules = sorted(summary.module.unique(), key=lambda m: int(m[1:]) if m.startswith('L') and m[1:].isdigit() else m)
    cohorts = sorted(summary.dataset.unique())
    fig, axes = plt.subplots(1, 2, figsize=(12, max(7, len(modules) * .105)), layout='constrained')
    for ax, column, title, cmap, low, high in [
        (axes[0], 'spearman_rho', 'PCA vs GSVA: within-cohort Spearman correlation', 'RdBu_r', -1, 1),
        (axes[1], 'coverage_fraction', 'Measured variable genes / module size', 'viridis', 0, 1)]:
        matrix = summary.pivot(index='module', columns='dataset', values=column).reindex(index=modules, columns=cohorts)
        im = ax.imshow(matrix, aspect='auto', interpolation='nearest', cmap=cmap, vmin=low, vmax=high)
        ax.set_xticks(range(len(cohorts)), cohorts, rotation=40, ha='right')
        labels = range(0, len(modules), 5) if len(modules) > 40 else range(len(modules))
        ax.set_yticks(list(labels), [modules[i] for i in labels], fontsize=8)
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, shrink=.4)
    fig.savefig(output / 'agreement_and_coverage.png', dpi=170)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 4), layout='constrained')
    for name in cohorts:
        values = summary.loc[summary.dataset.eq(name), 'spearman_rho'].dropna()
        ax.hist(values, bins=np.linspace(-1, 1, 26), histtype='step', label=f'{name} (n={len(values)})', linewidth=1.5)
    ax.set(xlabel='PCA versus GSVA Spearman correlation', ylabel='Modules', title='Agreement across samples within each cohort')
    ax.legend(fontsize=8)
    fig.savefig(output / 'agreement_distribution.png', dpi=170)
    plt.close(fig)
    examples = [m for m in ('L12', 'L19', 'L60', 'L166') if any(m in pair[0] for pair in scores.values())]
    if examples:
        fig, axes = plt.subplots(len(cohorts), len(examples), figsize=(3.2*len(examples), 2.4*len(cohorts)),
                                 squeeze=False, layout='constrained')
        for i, cohort in enumerate(cohorts):
            pca, gsva = scores[cohort]
            for j, module in enumerate(examples):
                ax = axes[i, j]
                if module in pca:
                    ax.scatter(pca[module], gsva[module], s=9, alpha=.5)
                ax.set(title=f'{cohort} / {module}', xlabel='PCA (SD)', ylabel='GSVA')
        fig.savefig(output / 'example_modules.png', dpi=170)
        plt.close(fig)


def run(datasets, modules_path, output, *, module_info=None, threads=1):
    modules = load_modules(modules_path)
    if not modules:
        raise ValueError('Empty module table')
    output.mkdir(parents=True, exist_ok=True)
    reports, summaries, all_scores, seen = [], [], {}, set()
    for processed, anndata in datasets:
        check = check_dataset(processed, anndata)
        if check['errors']:
            raise ValueError(f'{processed}: failed output checks: {check["errors"]}')
        expression, genes, clinical, metadata, transform = load_dataset(processed)
        accession = metadata.get('geo_series', metadata.get('study_accession', processed.name))
        if accession in seen:
            raise ValueError(f'Duplicate cohort {accession}: choose one snapshot per study')
        seen.add(accession)
        print(f'Scoring {processed.name}: {expression.shape[1]} samples, {expression.shape[0]} source features', flush=True)
        pca, gsva, summary, loadings, universe = score_expression(expression, genes, modules, threads=threads)
        summary.insert(0, 'dataset', processed.name)
        summary['accession'] = accession
        for kind, table in [('pca_scores', pca), ('gsva_scores', gsva)]:
            table.to_csv(output / f'{processed.name}_{kind}.csv', float_format='%.12g')
        loadings.to_csv(output / f'{processed.name}_pca_loadings.csv.gz', index=False, compression={'method': 'gzip', 'mtime': 0})
        measured = set(universe.index)
        pd.DataFrame([{'module': m, 'gene_symbol': g, 'measured_variable': g in measured}
                      for m, geneset in modules.items() for g in geneset]).to_csv(
                          output / f'{processed.name}_gene_coverage.csv.gz', index=False, compression={'method': 'gzip', 'mtime': 0})
        report = {'dataset': processed.name, 'accession': accession, 'samples': len(pca),
                  'variable_gene_background': len(universe), 'scored_modules': len(pca.columns),
                  'input_processed': str(processed.resolve()), 'input_anndata': str(anndata.resolve()),
                  'expression_transform': transform, 'output_check': check,
                  'input_sha256': {p.name: sha256(p) for p in processed.iterdir() if p.is_file()}}
        reports.append(report)
        all_scores[processed.name] = pca, gsva
        summaries.append(summary)
        write_json(output / f'{processed.name}_run.json', report)
    summary = pd.concat(summaries, ignore_index=True)
    if module_info is not None:
        info = pd.read_csv(module_info, dtype={'module': str})
        summary = summary.merge(info.drop(columns=['n_genes'], errors='ignore'), on='module', how='left', validate='many_to_one')
        summary['annotation_note'] = np.where(summary.module.eq('L60'), 'Source label retained; client questioned annotation; not resolved here', '')
    summary.to_csv(output / 'module_comparison.csv', index=False, float_format='%.12g')
    overview = summary.groupby('dataset').agg(modules=('module', 'size'), scored=('spearman_rho', 'count'),
        median_spearman=('spearman_rho', 'median'), median_coverage=('coverage_fraction', 'median'))
    overview['samples'] = [len(all_scores[name][0]) for name in overview.index]
    overview.to_csv(output / 'cohort_summary.csv', float_format='%.12g')
    figures(summary, output, all_scores)
    manifest = {'analysis': 'Within-cohort module PCA and GSVA comparison; no survival fits or pooled expression',
                'modules_sha256': sha256(modules_path), 'modules': len(modules), 'modules_path': str(modules_path.resolve()),
                'gsva': {'implementation': 'GSEApy', 'version': importlib.metadata.version('gseapy'),
                         'kcdf': 'Gaussian', 'weight': 1, 'mx_diff': True, 'abs_rnk': False, 'seed': 123,
                         'min_size': MIN_GENES, 'threads': threads},
                'gene_selection': 'One feature per exact symbol: highest mean log-expression, ties by feature ID; exclude constant genes; no imputation',
                'pca_orientation': 'PC1 sign follows mean standardized module expression; export orientation correlation and loadings',
                'versions': {p: importlib.metadata.version(p) for p in ['numpy', 'pandas', 'scipy', 'gseapy']},
                'code_sha256': code_hashes('ds_bu_bioinformatics.module_compare'), 'cohorts': reports}
    if module_info:
        manifest['module_info_sha256'] = sha256(module_info)
    write_json(output / 'analysis_manifest.json', manifest)
    lines = ['# PCA and GSVA: complete available cohorts', '',
             f'{len(reports)} cohorts; {sum(r["samples"] for r in reports)} expression profiles; {len(modules)} source modules.', '',
             'Scores are calculated separately within each cohort. No outcome associations, risk prediction, or pooling is performed.', '',
             '| Cohort | Profiles | Modules with both scores | Median Spearman | Median gene coverage |',
             '|---|---:|---:|---:|---:|']
    for name, row in overview.iterrows():
        lines.append(f'| {name} | {int(row.samples)} | {int(row.scored)}/{len(modules)} | {row.median_spearman:.3f} | {row.median_coverage:.1%} |')
    lines += ['', '## Interpretation and limits', '',
              '- Positive correlation means the methods order patients similarly within that cohort; it is not prognostic validation.',
              '- Modules with fewer than five measured variable genes are explicitly unscored, never filled with zero.',
              '- PCA is standardized within cohort; GSVA scores depend on the measured background and cohort composition. Neither is a calibrated clinical risk score.',
              '- Source gene symbols are matched exactly; no undocumented symbol alias mapping.',
              '- L60 retains the original workbook annotation; the label questioned by the client remains unresolved.',
              '- Taylor preprocessing is provisional: independent review found six modules whose PCA-GSVA correlation changed sign under alternate preprocessing. Primary-run reproducibility does not establish preprocessing robustness.',
              '- Complete structural exports establish computational eligibility, not uniform assay QC or independent biological validation.',
              '- This run includes available complete exports only. It does not reproduce the meeting claim of eight cohorts / 1,486 tumors.', '',
              '## Files', '', 'See module_comparison.csv for all module/cohort results, *_scores.csv for patient scores,',
              '*_pca_loadings.csv.gz for PC1 gene weights, *_gene_coverage.csv.gz for missing genes, and analysis_manifest.json for hashes and settings.', '',
              '![Agreement and gene coverage](agreement_and_coverage.png)', '',
              '![Correlation distributions](agreement_distribution.png)', '']
    (output / 'report.md').write_text('\n'.join(lines))
    (output / 'checksums.sha256').write_text(''.join(f'{sha256(p)}  {p.name}\n' for p in sorted(output.iterdir())
                                                 if p.is_file() and p.name != 'checksums.sha256'))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--modules', type=Path, required=True)
    parser.add_argument('--module-info', type=Path)
    parser.add_argument('--snapshot', action='append', required=True, help='Complete snapshot to include, repeatable; avoids historical duplicates')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error('threads must be positive')
    if any(not name or name in {'.', '..'} or Path(name).name != name for name in args.snapshot):
        parser.error('snapshot must be a single directory name')
    datasets = [(args.data_dir/'processed'/n, args.data_dir/'anndata'/n) for n in args.snapshot]
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error('Use a new empty output directory for each analysis run')
    if any(output.is_relative_to(p.resolve()) or p.resolve().is_relative_to(output) for pair in datasets for p in pair):
        parser.error('analysis output must not overlap dataset inputs')
    run(datasets, args.modules, args.output_dir, module_info=args.module_info, threads=args.threads)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
