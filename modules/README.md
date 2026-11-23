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
