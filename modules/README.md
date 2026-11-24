# Gene modules

`hallmark_placeholder.csv` holds **placeholder** modules for building and testing PCA module
scoring until the CampLab module gene lists are shared. Replace it (or pass another file with
`--modules`) when they are available. Format: one row per gene, with `gene` (HGNC symbol) and
`module` columns, the same format as `tcga-prad --module-file`.

| Module | MSigDB gene set | Genes |
|---|---|---:|
| `androgen_response` | HALLMARK_ANDROGEN_RESPONSE | 101 |
| `e2f_targets` | HALLMARK_E2F_TARGETS | 200 |
| `g2m_checkpoint` | HALLMARK_G2M_CHECKPOINT | 200 |
| `emt` | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 200 |
| `interferon_gamma` | HALLMARK_INTERFERON_GAMMA_RESPONSE | 200 |

Source: MSigDB Hallmark collection v2025.1.Hs, `h.all.v2025.1.Hs.symbols.gmt`
(https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt,
SHA-256 `f22066af72e215ccb7b89d88e492c07e1eef17534c2ca7b0f9902cfecbbdd8e9`, downloaded November 13, 2026).
Liberzon et al., *Cell Systems* 2015; MSigDB gene sets are distributed under CC BY 4.0.

E2F targets and G2M checkpoint are proliferation sets, which are known to be prognostic in
prostate cancer, so they serve as a positive control for the scoring method.

## Decipher gene set (`decipher_erho2013.csv`)

Comparison gene set for the Decipher genomic classifier (Veracyte), requested at the 11/16 client
meeting. Source: Erho et al., "Discovery and validation of a prostate cancer genomic classifier that
predicts early metastasis following radical prostatectomy", *PLoS One* 2013;8(6):e66855, Table 2
(https://pmc.ncbi.nlm.nih.gov/articles/PMC3691249/).

- The classifier has 22 markers. Table 2 names 19; rows 7, 10 and 17 have no name in the published
  table text and are not included.
- `marker_type` is the feature type from Table 2 (coding, 3' UTR, intronic, antisense, non-coding).
  Gene-level expression measures the whole gene, not these specific features.
- PCAT32 and GLYATL1P4 (PCAT-80) are non-coding and absent from the gene-level outputs. TSBP1 is
  stored as `C6orf10` in most cohorts and is matched through an alias in `decipher.py`.
- `erho_direction` (+1 / -1) is the sign of each gene's association with metastasis in Erho
  (GSE46691, the classifier's training cohort; 545 tumors, 212 with metastasis): Cohen's d of
  metastasis vs no metastasis, with a Welch p-value. Input: a gene-level summary (median of core
  probe sets per gene) of the GEO quantile-normalized deposit, after removing the genome-wide first
  principal component. That component explains 31% of variance, tracks overall sample intensity
  (FFPE quality), and is itself 0.28 SD higher in metastasis cases, so leaving it in biases the
  directions (LASP1, NFIB and PBX1 flip sign without the correction). EPPK1 and TSBP1 had no core
  probe set in that summary, so they have no direction. 10 of the 15 directions have p < 0.05.
  The summary was built outside the repository (November 16, 2026); recompute the directions from
  the `erho-gse46691` pipeline output with `decipher.directions_from_matrix` once that cohort is
  processed.
- The published paper gives neither per-marker directions nor model weights, so this file supports a
  gene-set comparison, not a reproduction of the Decipher score.

Do not pass this file to `module-scores`: `load_modules` drops rows with empty cells, and the
unsigned PC1 score is not meaningful for a gene set with mixed directions. Use `decipher-comparison`.
