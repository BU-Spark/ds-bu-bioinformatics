"""Decipher (Veracyte) gene-set comparison scores for each processed cohort.

The gene set is the named markers of the 22-marker genomic classifier in Erho et al., PLoS One
2013;8:e66855, Table 2 (modules/decipher_erho2013.csv). This is a gene-set comparison, NOT the
Decipher score: the commercial classifier is a locked random forest on exon-array features (many
intronic, UTR or antisense) whose weights are not published, three of its 22 markers are not named
in the source table, and two named markers are non-coding transcripts absent from gene-level data.

Two scores are written per cohort, both standardized within the cohort (mean 0, SD 1):
  * decipher_pca: the standard module score (PC1, higher = higher mean expression). Because the
    classifier mixes genes that rise and fall with metastasis, PC1 can follow either group, so its
    direction is not comparable between cohorts. Reported for completeness only.
  * decipher_signed: the mean of standardized expression after multiplying each gene by its
    direction of association with metastasis in Erho (GSE46691, the classifier's training cohort;
    stored in the gene-set file). Higher = more metastasis-like. Genes without a direction are left
    out of this score.

Cox models (per cohort, and pooled with study strata via pooled_analysis) compare these scores with
a reference module score, by default CampLab L12 (cycling kinase / proliferation).

Outputs (in --output-dir): <snapshot>_decipher_scores.csv, decipher_per_cohort.csv, decipher_pooled.csv.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from ds_bu_bioinformatics.module_scores import endpoint, load_dataset, module_matrix, pca_score, survival
from ds_bu_bioinformatics.pooled_analysis import RELIABLE, cohort_frame, fit_cox

# Gene-set symbol -> symbols used for the same gene in processed genes.csv files.
ALIASES = {'TSBP1': ('TSBP1', 'C6orf10', 'TSBP'), 'PCAT32': ('PCAT32', 'PCAT-32'),
           'GLYATL1P4': ('GLYATL1P4', 'PCAT80', 'PCAT-80')}
SCORES = ('decipher_pca', 'decipher_signed')


def load_gene_set(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    missing = {'gene', 'erho_direction'} - set(table.columns)
    if missing:
        raise ValueError(f'{path}: missing column(s) {sorted(missing)}')
    if not table.erho_direction.dropna().isin([-1, 1]).all():
        raise ValueError(f'{path}: erho_direction must be -1, 1 or empty')
    return table


def canonical_symbols(genes: pd.DataFrame) -> pd.DataFrame:
    """Copy of genes.csv with alias symbols renamed to the gene-set symbol."""
    genes = genes.copy()
    for symbol, names in ALIASES.items():
        genes.loc[genes.gene_symbol.isin(names), 'gene_symbol'] = symbol
    return genes


def directions_from_matrix(matrix: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    """Per-gene Cohen's d (cases minus controls), Welch p (normal approximation) and sign.

    Recomputes the erho_direction column from a genes x samples matrix and 0/1 metastasis labels."""
    labels = labels.reindex(matrix.columns)
    rows = []
    for gene, values in matrix.iterrows():
        a, b = values[labels == 1].astype(float), values[labels == 0].astype(float)
        d = (a.mean() - b.mean()) / math.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        t = (a.mean() - b.mean()) / math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
        rows.append({'gene': gene, 'erho_cohens_d': d, 'erho_welch_p': math.erfc(abs(t) / math.sqrt(2)),
                     'erho_direction': int(np.sign(d))})
    return pd.DataFrame(rows)


def signed_score(matrix: pd.DataFrame, directions: pd.Series) -> tuple[pd.Series, int]:
    """Mean of direction-signed standardized expression, rescaled to mean 0, SD 1."""
    used = [g for g in matrix.index if g in directions.index and not pd.isna(directions[g])]
    if not used:
        raise ValueError('no gene-set genes with a direction are measured')
    z = matrix.loc[used]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1, ddof=1), axis=0)
    score = z.mul(directions[used].astype(float), axis=0).mean(axis=0)
    return (score - score.mean()) / score.std(ddof=1), len(used)


def score_cohort(processed: Path, gene_set: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    expression, genes, clinical, metadata, transform = load_dataset(processed)
    matrix = module_matrix(expression, canonical_symbols(genes), gene_set.gene)
    pca, variance, _ = pca_score(matrix)
    signed, used = signed_score(matrix, gene_set.set_index('gene').erho_direction)
    scores = pd.DataFrame({'decipher_pca': pca, 'decipher_signed': signed}).reindex(clinical.sample_id)
    scores.index.name = 'sample_id'
    return scores, {'dataset': processed.name, 'genes_pca': len(matrix), 'genes_signed': used,
                    'pc1_variance': round(variance, 4), 'expression_transform': transform}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--data-dir', type=Path, default=Path('data'), help='Folder containing processed/<snapshot>')
    parser.add_argument('--gene-set', type=Path, default=Path('modules/decipher_erho2013.csv'))
    parser.add_argument('--scores-dir', type=Path, help='module-scores output with the reference module; '
                        'default <data-dir>/analysis/module_scores')
    parser.add_argument('--reference-module', default='L12', help='Module compared with the Decipher gene set')
    parser.add_argument('--dataset', action='append', default=[], help='Snapshot to score (repeatable); default all')
    parser.add_argument('--output-dir', type=Path, help='Default: <data-dir>/analysis/decipher')
    args = parser.parse_args(argv)
    gene_set = load_gene_set(args.gene_set)
    scores_dir = args.scores_dir or args.data_dir / 'analysis' / 'module_scores'
    output = args.output_dir or args.data_dir / 'analysis' / 'decipher'
    output.mkdir(parents=True, exist_ok=True)
    processed_root = args.data_dir / 'processed'
    datasets = args.dataset or sorted(p.name for p in processed_root.iterdir() if (p / 'clinical.csv').exists())

    per_cohort, frames = [], []
    for name in datasets:
        processed = processed_root / name
        scores, info = score_cohort(processed, gene_set)
        scores.to_csv(output / f'{name}_decipher_scores.csv')
        reference = scores_dir / f'{name}_scores.csv'
        if reference.exists():
            module_scores = pd.read_csv(reference, dtype={'sample_id': str}).set_index('sample_id')
            if args.reference_module in module_scores:
                scores[args.reference_module] = module_scores[args.reference_module].reindex(scores.index)
        clinical = pd.read_csv(processed / 'clinical.csv', dtype={'sample_id': str}, low_memory=False)
        metadata = json.loads((processed / 'metadata.json').read_text())
        try:
            endpoint_name, event, months, eligible = endpoint(clinical, metadata, name)
        except ValueError as error:
            per_cohort.append({**info, 'note': str(error)})
            continue
        info.update({'endpoint': endpoint_name, 'eligible': int(eligible.sum()), 'events': int(event[eligible].eq(1).sum()),
                     'corr_pca_vs_signed': float(scores.decipher_pca.corr(scores.decipher_signed))})
        if args.reference_module in scores:
            info[f'corr_signed_vs_{args.reference_module}'] = float(scores.decipher_signed.corr(scores[args.reference_module]))
        for column in SCORES:
            result = survival(scores[column][eligible], months[eligible], event[eligible])
            info.update({f'{column}_{k}': v for k, v in result.items()})
        per_cohort.append(info)
        if name.startswith(RELIABLE) and args.reference_module in scores:
            frames.append(cohort_frame(processed, scores))
    per = pd.DataFrame(per_cohort)
    per.to_csv(output / 'decipher_per_cohort.csv', index=False)

    pooled_rows = []
    if frames:
        pooled = pd.concat(frames)
        ref = args.reference_module
        adjustments = {'unadjusted': [], 'adj_gleason': ['gleason_cat'], 'adj_gleason_pt3': ['gleason_cat', 'pt3_plus']}
        for terms in (['decipher_signed'], ['decipher_pca'], [ref], ['decipher_signed', ref]):
            for model, covariates in adjustments.items():
                pooled_rows.append({'terms': ' + '.join(terms), 'model': model,
                                    **fit_cox(pooled, terms, covariates, 'study')})
    pool = pd.DataFrame(pooled_rows)
    pool.to_csv(output / 'decipher_pooled.csv', index=False)

    with pd.option_context('display.width', 200, 'display.max_columns', 30):
        columns = [c for c in ['dataset', 'endpoint', 'genes_signed', 'events', 'decipher_signed_hazard_ratio_per_sd',
                               'decipher_signed_hr_ci_low', 'decipher_signed_hr_ci_high', 'decipher_signed_p_value',
                               'decipher_pca_hazard_ratio_per_sd'] if c in per]
        print(per[columns].round(3).to_string(index=False))
        if not pool.empty:
            show = ['terms', 'model', 'n', 'events', 'cohorts'] + [c for c in pool if c.endswith(('_hr', '_p'))]
            print('\nPooled, stratified by study:')
            print(pool[show].round(4).to_string(index=False))
    print(f'\nWrote {output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
