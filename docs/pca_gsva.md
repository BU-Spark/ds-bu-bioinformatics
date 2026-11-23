# PCA and GSVA on complete exports

This workflow compares module scores across patients **within each cohort**, using the
actual CampLab Supplementary Table 2 modules. It does not fit Cox models, pool expression
across platforms, change L60's annotation, or reproduce the clinical Decipher algorithm.

```bash
uv sync --locked
uv run python -m ds_bu_bioinformatics.module_compare \
  --modules modules/camplab_modules.csv \
  --module-info modules/camplab_modules_info.csv \
  --snapshot taylor-v2 --snapshot cambridge-v2 --snapshot stockholm-v1 \
  --snapshot long-v1 --snapshot mortensen-v1 \
  --output-dir data/analysis/camplab-pca-gsva-v1 --threads 2
```

Use a new empty output directory for every run. Each selected snapshot must pass the
full CSV/H5AD structural checker. Select one snapshot per study; duplicate cohorts are
rejected. Add another `--snapshot` only after its exports are complete. Incomplete
acquisitions and raw-only caches are excluded, not silently scored or counted.

## Shared preprocessing

1. Preserve exported source values. For scoring only, transform raw counts to
   `log2(CPM+1)` using frozen full-source library totals; transform TPM/FPKM and declared
   linear microarrays with `log2(x+1)`. Deposited log microarrays remain unchanged.
   Legacy exports without a declared kind retain the documented maximum-value heuristic;
   the exact applied transform is recorded in each run JSON. Review this legacy heuristic
   before using another undeclared source scale.
2. Use exact source gene symbols. Choose one feature per symbol by highest mean
   transformed expression, breaking ties by feature ID. This uses no outcome data.
   Remove constant genes. Both methods use the same selected genes; GSVA uses the
   **whole variable-gene background**, not just module genes.
3. PCA standardizes genes within cohort and scores PC1, oriented to correlate positively
   with mean standardized module expression. Scores have mean zero and sample SD one.
   Export variance explained, loadings, orientation correlation and positive-loading share.
   Constant mean gene-z makes orientation undefined and is reported explicitly.
4. GSVA uses pinned **GSEApy 1.3.0**, Gaussian kernel on continuous log-expression,
   weight=1, mx_diff=True, abs_rnk=False, seed=123. This is GSVA, not ssGSEA.
   Both methods require at least five measured variable module genes. Smaller modules
   remain explicitly unscored; all 166 modules remain in the coverage table.

The five-gene threshold is an exploratory choice, not evidence of reliable scores for
poorly covered modules. Examine coverage, PC1 variance and orientation before interpretation.
GSEApy is a Python implementation; this run does not claim bitwise equivalence to an R
Bioconductor GSVA version.

## Outputs

- `report.md`, cohort summary, agreement/coverage figures and example module scatterplots.
- Per-cohort PCA and GSVA patient scores, PC1 gene loadings, and per-gene coverage audits.
- `module_comparison.csv`: every module × cohort, missing-score reasons, gene coverage,
  PC1 QC, Pearson and Spearman agreement. These correlations are method comparisons,
  not associations with outcome or evidence of independent prognostic value.
- `analysis_manifest.json`: input hashes, implementation versions, source paths, code
  hashes, full structural checks and settings; `checksums.sha256` covers results.

Source workbook identity is frozen in `modules/camplab_modules_source.json`. The original
annotation for L60 remains unchanged and explicitly flagged for review. Hallmark controls
are separate from the 166-module analysis.

Method references: [GSVA vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html)
and [GSEApy GSVA API](https://gseapy.readthedocs.io/en/master/_modules/gseapy.html).

## Taylor sensitivity caveat

Independent review found six Taylor modules whose PCA–GSVA correlation changed sign
under alternate preprocessing. The legacy Taylor transformation is provisional;
primary-run reproducibility does not demonstrate preprocessing robustness. Retain
this caveat when sharing the tables or interpreting individual correlation directions.
