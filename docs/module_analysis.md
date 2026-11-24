# Module analysis: FDR, clinical adjustment and cross-cohort meta-analysis

`module-analysis` adds statistics on top of the PCA module scores written by
[`module-scores`](module_scoring.md). It does not change any pipeline output or
`module_scores.py`; it reads their outputs and the processed clinical tables.

Code: `src/ds_bu_bioinformatics/module_analysis.py`. Tests: `tests/test_module_analysis.py`
(synthetic data only, no downloads).

## How to run

The CampLab modules come from Supplementary Table 2 of the single-cell atlas
(`Supp_Table_2_Module_List.xlsx`). Convert them once, score them, then run the analysis:

```bash
python -m ds_bu_bioinformatics.camplab_modules --input path/to/Supp_Table_2_Module_List.xlsx
python -m ds_bu_bioinformatics.module_scores --modules modules/camplab_modules.csv --output-dir data/analysis/camplab_scores
python -m ds_bu_bioinformatics.module_analysis --scores-dir data/analysis/camplab_scores --module-info modules/camplab_modules_info.csv
```

`camplab_modules` writes `modules/camplab_modules.csv` (`gene,module`, the `module-scores` input)
and `modules/camplab_modules_info.csv` (`module,category,annotation,n_genes`). It parses the
workbook with the standard library, checks that the header row holds module IDs (`L1`, `L2`, ...),
keeps gene symbols as written and reports anything it removes. The source table has 166 modules and
12,680 genes; each gene belongs to one module. 28 modules have fewer than 5 genes and are skipped by
`module-scores` (`MIN_GENES = 5`).

Options for `module_analysis`:

| Option | Default | Meaning |
|---|---|---|
| `--data-dir` | `data` | Root with `processed/<snapshot>/` |
| `--scores-dir` | `data/analysis/module_scores` | `module-scores` output (`summary.csv`, `<snapshot>_scores.csv`) |
| `--module-info` | none | Labels such as "L12 · Cycling kinase" in tables and figures |
| `--exclude` | none | Leave a snapshot out of the meta-analysis (repeatable), e.g. `--exclude taylor-v2` |
| `--top` | 20 | Modules shown in the report and heatmap |
| `--output-dir` | `<scores-dir>/analysis` | Where results are written |

On Windows machines where Smart App Control blocks the installed `.exe` launchers, use the
`python -m ds_bu_bioinformatics.<module>` form shown above.

## Outputs

| File | Contents |
|---|---|
| `per_dataset.csv` | One row per dataset × module: unadjusted results from `summary.csv` plus `q_value`, adjusted HR/CI/p, `adj_q_value`, `added_value_lr_p`, `delta_c_index` |
| `clinical_models.csv` | Per dataset: endpoint, covariates used and dropped, complete cases, events, events per parameter, clinical-only C-index |
| `meta.csv` | One row per module: pooled HR, CI, p, q, I², τ², number of datasets and how many agree in direction, for unadjusted and adjusted (`adj_`) estimates |
| `analysis.md` | Readable summary of the three tables above |
| `heatmap.png` | Unadjusted HR per SD for the top modules × datasets; blue = lower risk, red = higher risk, `*` = q < 0.05 |
| `forest_top.png` | Per-dataset HR and 95% CI with the pooled estimate for the top 8 modules |

Outputs are under `data/`, which Git ignores.

## Method

1. **Within-dataset FDR.** Unadjusted results are taken unchanged from `module-scores`
   (`summary.csv`). Benjamini-Hochberg q-values are computed within each dataset across its scored
   modules.
2. **Clinical-adjusted Cox.** For each dataset and module, a Cox model of the dataset's own primary
   endpoint (the same `endpoint()` and eligibility used by `module-scores`) on the module score plus
   clinical covariates, using complete cases among eligible records. The score's added value is
   measured against a clinical-only model fitted on the same patients: a likelihood-ratio test
   (1 df) and the change in C-index. BH q-values are computed within each dataset.
3. **Random-effects meta-analysis.** For each module, the log hazard ratios per SD from all datasets
   are pooled with DerSimonian-Laird, separately for unadjusted and adjusted estimates. Standard
   errors are derived from the 95% CIs. Reported: pooled HR, 95% CI, p, BH q across modules, I²,
   τ², and the number of datasets whose estimate has the same sign as the pooled one.

Endpoints are not relabeled: TCGA uses PFI, Taylor DFS, Cambridge and Stockholm BCR, and Belfast
metastasis. The pooled estimate therefore answers "does the module track disease progression in
general", not one endpoint.

### Clinical covariates

Each cohort names and codes covariates differently, so the mapping is explicit in `COVARIATES`:

| Covariate | Coding | TCGA-PRAD | Taylor | Cambridge | Stockholm | Belfast |
|---|---|---|---|---|---|---|
| `grade_group` | ISUP grade group 1–5 from Gleason primary + secondary (3+4 = 2, 4+3 = 3) | `legacy_patient__GLEASON_PATTERN_PRIMARY/SECONDARY` | `sample__GLEASON_SCORE_1/2` | `geo__tumour_gleason` (`7=3+4`) | `geo__tumour_gleason` | — |
| `gleason_8_plus` | 1 for total Gleason 8–10, 0 for ≤7; used where only the total score is recorded | — | — | — | — | `geo__gleason_grade` |
| `pt3_plus` | 1 for pathological T3/T4, 0 for T0–T2, missing for Tx | `path_t_group` | `patient__PATH_T_STAGE` | `geo__pathology_stage` | `geo__pathology_stage` | — |
| `ct3_plus` | 1 for clinical T3/T4, 0 for T1–T2, missing for Unknown | — | — | — | — | `geo__t_stage` |
| `age` | Years | `age` | — | `geo__age_at_diag` | — | `geo__patient_age_years` |
| `n1` | 1 for N1, 0 for N0, missing for Nx | `patient__PATH_N_STAGE` | — | — | — | — |

Belfast is a biopsy cohort, so it records only the total Gleason score (no grade group) and clinical
rather than pathological T stage; see [belfast.md](belfast.md).

Rules:

- Missing-value tokens (`unknown`, `NA`, blank, ...) become missing. Any other unrecognized value
  stops the run instead of silently becoming missing; for example, Stockholm's `T0MxNx` was found
  this way and is coded as not pT3+ (no residual tumour in the specimen).
- A covariate is dropped for a dataset if more than 20% of eligible records lack it, or if it has
  only one value. Dropped covariates are listed in `clinical_models.csv`.
- Fewer than 10 events per model parameter is flagged as possibly unstable; the model is still fit.
- Cohorts not listed in `COVARIATES` (currently Sboner, Mortensen, Long and CIT) get unadjusted
  results and are included in the unadjusted meta-analysis only. Adding a cohort means adding its
  verified column mapping.

## First results with the CampLab modules

Run on November 16, 2026 on `tcga-prad-v2`, `taylor-v2`, `cambridge-v2`, `stockholm-v1` and
the earlier Belfast metastasis export (130–135 of 166 modules scored per dataset).
The current Belfast pipeline uses different feature selection and defaults to BCR.
These historical estimates have not been reproduced with the current export.

Clinical models:

| Dataset | Endpoint | Covariates | Complete cases (events) | Events/parameter | Clinical C-index |
|---|---|---|---:|---:|---:|
| tcga-prad-v2 | PFI | age, grade group, pT3+, N1 | 409 (82) | 16.4 | 0.72 |
| taylor-v2 | DFS | grade group, pT3+ | 130 (27) | 9.0 ⚠ | 0.83 |
| cambridge-v2 | BCR | age, grade group, pT3+ | 111 (19) | 4.8 ⚠ | 0.76 |
| stockholm-v1 | BCR | grade group, pT3+ | 89 (44) | 14.7 | 0.76 |
| belfast-v1 | Metastasis | age, Gleason 8+, clinical T3+ | 223 (22) | 5.5 ⚠ | 0.61 |

Modules with BH q < 0.05:

| Dataset | Unadjusted | Adjusted |
|---|---:|---:|
| tcga-prad-v2 | 29 | 0 |
| stockholm-v1 | 21 | 0 |
| belfast-v1 | 2 | 0 |
| cambridge-v2 | 1 | 0 |
| taylor-v2 | 0 | 0 |

Top modules by random-effects meta-analysis (all 5 datasets):

| Module | Pooled HR per SD (95% CI) | q | I² | Same direction | Adjusted pooled HR (95% CI) | Adjusted q |
|---|---|---:|---:|---:|---|---:|
| L60 · Neuroendocrine | 0.56 (0.49–0.64) | 2.5e-14 | 0% | 5/5 | 0.75 (0.62–0.92) | 0.30 |
| L12 · Cycling kinase | 1.61 (1.42–1.82) | 4.3e-12 | 10% | 5/5 | 1.35 (1.13–1.59) | 0.087 |
| L11 · DNA replication | 1.52 (1.34–1.74) | 1.0e-08 | 9% | 5/5 | 1.22 (1.06–1.41) | 0.30 |
| L16 · Metallothionein | 0.67 (0.59–0.77) | 9.6e-07 | 0% | 5/5 | 0.86 (0.74–0.99) | 0.59 |
| L59 · AR pathway | 0.70 (0.62–0.80) | 6.1e-06 | 0% | 5/5 | 0.86 (0.74–0.99) | 0.59 |
| L72 · Luminal | 0.72 (0.63–0.82) | 3.3e-05 | 0% | 5/5 | 0.88 (0.76–1.02) | 0.73 |
| L63 · Stromal remodelling | 0.74 (0.65–0.84) | 4.1e-05 | 0% | 5/5 | 0.89 (0.77–1.02) | 0.73 |
| L64 · AR pathway | 0.74 (0.65–0.84) | 0.00015 | 0% | 5/5 | 0.88 (0.76–1.02) | 0.75 |
| L61 · PI3K pathway | 0.75 (0.65–0.85) | 0.00032 | 0% | 4/5 | 0.86 (0.74–1.00) | 0.59 |
| L58 · AR pathway | 0.66 (0.54–0.80) | 0.00037 | 38% | 5/5 | 0.84 (0.64–1.08) | 0.83 |

The full top-20 table, the heatmap and the forest plots are in the run's `analysis/` folder.

### Findings

- **Proliferation predicts worse outcome in every cohort.** L12 (cycling kinase) has a pooled HR of
  1.61 per SD across 5/5 cohorts and is the strongest module after adjusting for Gleason and stage
  (adjusted HR 1.35, 95% CI 1.13–1.59, adjusted q = 0.087). L11 (DNA replication) behaves similarly.
  This matches the Hallmark E2F/G2M positive control in [module_scoring.md](module_scoring.md).
- **AR pathway, luminal and metallothionein modules are associated with lower risk** in 5/5 cohorts
  (pooled HR about 0.70).
- **L60 ("Neuroendocrine") is protective in all five cohorts** (pooled HR 0.56, I² = 0%), including
  Belfast's metastasis endpoint. Neuroendocrine features are usually linked to aggressive disease, so
  the module's genes and label should be reviewed with CampLab.
- **Most of the signal overlaps with Gleason grade and stage.** After adjustment, no module is
  significant at q < 0.05 within a single cohort.
- **Belfast (metastasis)** is the exception for proliferation: L12 and L11 do not predict metastasis
  there, while L58 (AR pathway) and L60 do, and add to the clinical model (C-index 0.61 → 0.76 for
  L58). See [belfast.md](belfast.md).

### Caveats

- Cambridge (19 events), Taylor (27) and Belfast (22) have fewer than 10 events per parameter in the
  adjusted models; their adjusted estimates are unstable.
- Taylor's expression scale is still unresolved (see [module_scoring.md](module_scoring.md)); rerun
  with `--exclude taylor-v2` to check that conclusions do not depend on it.
- Endpoints differ by cohort (PFI, DFS, BCR, metastasis), so pooled estimates mix related but
  different outcomes.
- Belfast covariates are coarser: total Gleason only and clinical rather than pathological T stage.
- C-indices are in-sample, not held-out performance.
- PCA scoring uses one probe per gene (highest mean), as implemented in `module-scores`.

## Open questions

1. Can `modules/camplab_modules.csv` be committed, given that the atlas paper is unpublished?
   Until confirmed it stays local; anyone with the workbook can regenerate it with `camplab_modules`.
2. Should L60 and the other "Neuroendocrine" modules be reviewed with CampLab?
3. Add covariate mappings for Sboner, Mortensen, Long and CIT once their columns are verified.
4. Next step: biomarker models (multi-module combinations vs. single genes, existing signatures and
   clinical variables) with held-out validation across cohorts.