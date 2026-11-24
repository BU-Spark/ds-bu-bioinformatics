"""PCA module scores for every processed cohort, including cohorts without usable outcomes.

module-scores computes PC1 scores and immediately fits a Cox model, so it cannot score cohorts with no
eligible outcomes. This command computes the scores only, using module-scores' own functions
(module_matrix, pca_score, MIN_GENES), so scores are identical to module-scores for the same input.
It does not modify module_scores.py or any processed export.

Preparation before PCA follows each cohort's declared expression kind (metadata.json):
  * raw_counts  -> log2(CPM + 1) using the full-source library sizes stored by the pipeline
  * linear, RPKM/FPKM, or any matrix with max > 50 -> log2(x + 1) (the module-scores rule)
  * log-scale microarray / RNA-seq -> used as deposited
Only feature rows of module genes are read from large expression files (chunked), so exon-array cohorts fit
in memory. Outputs: <snapshot>_scores.csv, coverage.csv, cohorts.csv and summary.md.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.module_scores import LOG_SCALE_MAX, MIN_GENES, load_modules, module_matrix, pca_score

LINEAR_KINDS = {'linear_microarray', 'rpkm', 'fpkm', 'tpm'}
CHUNK_ROWS = 50_000


def read_module_rows(path: Path, wanted: set[str]):
    """Return (rows for wanted feature IDs, max value over the whole matrix, total feature rows)."""
    parts, maximum, total = [], -np.inf, 0
    for chunk in pd.read_csv(path, index_col=0, dtype={'feature_id': str}, chunksize=CHUNK_ROWS):
        chunk.index = chunk.index.astype(str)
        total += len(chunk)
        maximum = max(maximum, float(np.nanmax(chunk.to_numpy())))
        keep = chunk.loc[chunk.index.isin(wanted)]
        if len(keep):
            parts.append(keep)
    rows = pd.concat(parts) if parts else pd.DataFrame()
    return rows, maximum, total


def prepare(expression: pd.DataFrame, maximum: float, metadata: dict, dataset: str):
    """Apply the declared, documented transform; return (matrix, description)."""
    info = metadata.get('expression', {}) or {}
    kind = info.get('kind')
    if kind == 'raw_counts':
        sizes = info.get('library_sizes')
        if not isinstance(sizes, dict):
            raise ValueError(f'{dataset}: raw counts without stored library sizes')
        missing = [s for s in expression.columns if s not in sizes]
        if missing:
            raise ValueError(f'{dataset}: library sizes missing for {missing[:3]}')
        totals = pd.Series({s: float(sizes[s]) for s in expression.columns})
        if (totals <= 0).any():
            raise ValueError(f'{dataset}: non-positive library size')
        return np.log2(expression.clip(lower=0).div(totals, axis=1) * 1e6 + 1), 'log2(CPM + 1) from full-source library sizes'
    if kind in LINEAR_KINDS or maximum > LOG_SCALE_MAX:
        return np.log2(expression.clip(lower=0) + 1), f'log2(x + 1) ({kind or "max > " + str(int(LOG_SCALE_MAX))})'
    return expression, 'none (already log scale)'


def eligible_outcomes(clinical: pd.DataFrame) -> str:
    found = []
    for column in clinical.columns:
        if column.endswith('_eligible') and '__' not in column and f'{column[:-9]}_event' in clinical:
            n = int(clinical[column].astype(str).str.lower().eq('true').sum())
            if n:
                found.append(f'{column[:-9]} ({n})')
    return ', '.join(found)


def score_cohort(processed: Path, modules: dict[str, list[str]]):
    name = processed.name
    metadata = json.loads((processed / 'metadata.json').read_text())
    genes = pd.read_csv(processed / 'genes.csv', dtype=str, keep_default_na=False)
    clinical = pd.read_csv(processed / 'clinical.csv', dtype={'sample_id': str}, low_memory=False)
    if 'gene_symbol' not in genes:
        raise ValueError(f'{name}: genes.csv has no gene_symbol column')
    every = {g for symbols in modules.values() for g in symbols}
    wanted = set(genes.loc[genes.gene_symbol.isin(every), 'feature_id'])
    rows, maximum, total = read_module_rows(processed / 'expression.csv.gz', wanted)
    if total != len(genes):
        raise ValueError(f'{name}: expression and genes.csv row counts differ; run check-outputs')
    if rows.empty:
        raise ValueError(f'{name}: no module genes measured')
    if rows.columns.tolist() != clinical.sample_id.tolist():
        raise ValueError(f'{name}: expression columns and clinical.csv samples differ; run check-outputs')
    expression, transform = prepare(rows, maximum, metadata, name)
    sub_genes = genes.set_index('feature_id').loc[expression.index].reset_index()
    scores, coverage = {}, []
    for module, symbols in modules.items():
        matrix = module_matrix(expression, sub_genes, symbols)
        record = {'dataset': name, 'module': module, 'module_genes': len(symbols), 'genes_used': len(matrix)}
        if len(matrix) < MIN_GENES:
            coverage.append({**record, 'scored': False, 'note': f'fewer than {MIN_GENES} module genes measured'})
            continue
        score, variance, _ = pca_score(matrix)
        scores[module] = score
        coverage.append({**record, 'scored': True, 'pc1_variance_explained': round(variance, 4)})
    table = pd.DataFrame(scores, index=expression.columns)
    table.index.name = 'sample_id'
    cov = pd.DataFrame(coverage)
    scored = cov[cov.scored]
    cohort = {'dataset': name, 'samples': expression.shape[1], 'feature_rows': total,
              'unique_genes': int(genes.gene_symbol.replace('', np.nan).nunique()), 'transform': transform,
              'modules_scored': len(scored),
              'median_pc1_variance': round(float(scored.pc1_variance_explained.median()), 4) if len(scored) else np.nan,
              'eligible_outcomes': eligible_outcomes(clinical) or 'none'}
    return table, cov, cohort


def markdown(cohorts: pd.DataFrame, failures: list[tuple[str, str]], n_modules: int) -> str:
    lines = ['# PCA module scores, all cohorts', '',
             f'{len(cohorts)} cohorts scored with {n_modules} modules (PC1 per module, module-scores method).', '',
             '| Cohort | Samples | Genes | Transform before PCA | Modules scored | Median PC1 variance | Outcomes (eligible) |',
             '|---|---:|---:|---|---:|---:|---|']
    for r in cohorts.itertuples():
        var = '' if pd.isna(r.median_pc1_variance) else f'{r.median_pc1_variance:.0%}'
        lines.append(f'| {r.dataset} | {r.samples} | {r.unique_genes:,} | {r.transform} | {r.modules_scored} | {var} | '
                     f'{r.eligible_outcomes} |')
    lines.append(f'| **Total** | **{int(cohorts.samples.sum()):,}** | | | | | |')
    if failures:
        lines += ['', '## Not scored', ''] + [f'- {name}: {error}' for name, error in failures]
    lines += ['', 'Scores are scaled to mean 0, SD 1 within each cohort; higher = higher module expression.',
              'Outcome testing is separate (module-scores / module_analysis) and applies only to cohorts with eligible outcomes.']
    return '\n'.join(lines) + '\n'


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--data-dir', type=Path, default=Path('data'))
    parser.add_argument('--modules', type=Path, default=Path('modules/camplab_modules.csv'))
    parser.add_argument('--snapshot', action='append', help='Score only this snapshot (repeatable); default: all')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    modules = load_modules(args.modules)
    root = args.data_dir / 'processed'
    names = args.snapshot or sorted(p.name for p in root.iterdir()
                                    if p.is_dir() and (p / 'validation.json').is_file())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cohorts, coverage, failures = [], [], []
    for name in names:
        try:
            table, cov, cohort = score_cohort(root / name, modules)
        except Exception as error:   # report and continue with the remaining cohorts
            failures.append((name, str(error)))
            print(f'{name}: NOT SCORED - {error}')
            continue
        table.to_csv(args.output_dir / f'{name}_scores.csv', float_format='%.6g')
        coverage.append(cov)
        cohorts.append(cohort)
        print(f"{name}: {cohort['samples']} samples, {cohort['modules_scored']} modules scored ({cohort['transform']})")
    if not cohorts:
        print('No cohort was scored.')
        return 1
    cohorts = pd.DataFrame(cohorts)
    pd.concat(coverage, ignore_index=True).to_csv(args.output_dir / 'coverage.csv', index=False)
    cohorts.to_csv(args.output_dir / 'cohorts.csv', index=False)
    report = markdown(cohorts, failures, len(modules))
    (args.output_dir / 'summary.md').write_text(report, encoding='utf-8')
    print('\n' + report)
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
