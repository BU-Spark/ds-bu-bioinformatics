# Pooled Cox with a study term, and the Decipher gene-set comparison

Follow-up to the 11/16 client meeting. Joshua asked for (1) a combined-cohort Cox model with a study
covariate, (2) the PC1 orientation rule to be verified and documented (see
[module_scoring.md](module_scoring.md#pc1-orientation-rule)), and (3) the Decipher gene set (Veracyte)
as a comparison signature. Results below are from November 24, 2026, on the CampLab module scores
and the processed snapshots of November 15.

## Commands

```bash
# Pooled Cox for every module, plus the PC1 orientation check
uv run pooled-analysis --data-dir data --scores-dir data/analysis/module_scores \
    --modules <camplab gene,module CSV> --module-info <module,category,annotation CSV>

# Decipher gene-set scores per cohort, per-cohort Cox, pooled comparison with L12
uv run decipher-comparison --data-dir data --scores-dir data/analysis/module_scores
```

`--data-dir` is the folder that contains `processed/<snapshot>`. Outputs go to
`<scores-dir>/pooled/` and `<data-dir>/analysis/decipher/`.

## Pooled model with a study term

All eligible patients of the five reliable cohorts (TCGA, Stockholm, Mortensen, Sboner, Cambridge)
in one Cox model: **1,005 patients, 342 events**. Module scores are already standardized within each
cohort, so hazard ratios are per SD of score. Each cohort keeps its own endpoint.

| Model | Cohorts | Patients / events | Modules with q < 0.05 (of 136) |
|---|---:|---|---|
| Stratified by study (primary) | 5 | 1,005 / 342 | **54** (36 higher risk, 18 lower risk) |
| Study as dummy covariates (Joshua's request) | 5 | 1,005 / 342 | 54 |
| Stratified, adjusted for Gleason (≤6 / 7 / 8+) | 5 | 1,003 / 342 | **15** |
| Stratified, adjusted for Gleason + pT3+ | 4 (Sboner has no stage) | 713 / 175 | **2**: L12, L11 |

- **The study term changes little.** Stratifying by study and adding study dummies agree. Both are
  close to the fixed-effect summary presented on 11/16 (48). `module-analysis` uses a random-effects
  meta-analysis, which allows true differences between cohorts and gives 25. The gap is the modules
  whose effect varies between cohorts.
- **Some modules do survive clinical adjustment** once patients are pooled. This differs from the
  earlier finding that nothing survived, which came from Belfast alone and from a meta-analysis with
  only three adjusted cohorts.
  - L12 (cycling kinase): HR 1.40 (1.21–1.62) after Gleason + pT3+; 1.30 (1.18–1.44) after Gleason.
  - L11 (DNA replication): HR 1.31 after Gleason + pT3+. L11's PC1 is weakly coherent in one cohort
    (correlation 0.29 with mean module expression), so interpret with care.
  - After Gleason alone: L11, L12, L129, L140, L144, L148, L149, L151, L157 and L158 (higher risk);
    L58, L59, L60, L61 and L72 (AR, luminal and PI3K, lower risk).
- Gleason is coarsened to three levels so cohorts that record only the total score (Mortensen,
  Sboner) can be included. The stage-adjusted model loses Sboner, which carries 165 of the 342 events.

## Decipher gene-set comparison

Gene set and direction source: [modules/README.md](../modules/README.md#decipher-gene-set-decipher_erho2013csv).
17 of the 19 named markers are measured in every cohort except Sboner (8). The signed score uses the
15 with a direction learned in Erho (Decipher's training cohort, not used here for testing).

**This is a gene-set comparison, not the Decipher score.** Decipher is a locked random forest on
exon-array features with unpublished weights, three of its 22 markers are unnamed in the source
table, and two named markers are non-coding.

The standard PCA score is not usable for this gene set. The genes go in both directions with risk,
so PC1 follows either group and points the wrong way in Cambridge (HR 0.51) and Taylor (0.61).

### Signed score per cohort

| Cohort | Endpoint | Events | HR per SD (95% CI) | p |
|---|---|---:|---|---:|
| TCGA | PFI | 91 | 1.80 (1.49–2.16) | <0.001 |
| Stockholm | BCR | 45 | 1.69 (1.27–2.26) | <0.001 |
| Mortensen | BCR | 22 | 1.92 (1.25–2.95) | 0.003 |
| Sboner (8 genes) | Prostate-cancer death | 165 | 1.16 (1.01–1.34) | 0.035 |
| Cambridge | BCR | 19 | 1.84 (1.10–3.08) | 0.020 |
| Taylor (provisional) | DFS | 27 | 2.11 (1.45–3.07) | <0.001 |
| Long (provisional) | BCR | 40 | 1.35 (1.06–1.72) | 0.015 |
| Belfast | Metastasis | 22 | 1.58 (1.09–2.28) | 0.016 |

Higher Decipher-like expression goes with worse outcome in all eight cohorts. That includes Long,
where no CampLab module has shown a signal, and Belfast's metastasis endpoint, which is what
Decipher was built to predict.

### Pooled with L12 (stratified by study, five reliable cohorts)

| Model | Decipher gene set HR | L12 HR | C-index |
|---|---|---|---:|
| Decipher alone, unadjusted | 1.42 (1.29–1.56) | | 0.622 |
| L12 alone, unadjusted | | 1.55 (1.42–1.69) | 0.639 |
| Decipher alone + Gleason + pT3+ | 1.36 (1.16–1.59) | | 0.728 |
| L12 alone + Gleason + pT3+ | | 1.40 (1.21–1.62) | 0.730 |
| **Both + Gleason + pT3+** | **1.21 (1.01–1.44), p = 0.04** | **1.29 (1.10–1.53), p = 0.002** | 0.735 |

Both stay significant in the same adjusted model, so the Decipher gene set and L12 carry partly
separate information. They correlate 0.21–0.67 across cohorts. C-indices are in-sample and
stratified, so they are not prediction estimates. A held-out comparison is part of the
prediction-model step.

## Caveats

- The endpoints differ (PFI, BCR, cancer death). Stratification gives each cohort its own baseline
  hazard but still assumes a common hazard ratio across endpoints.
- These are first-pass results on the November 15 snapshots. Taylor and Long remain provisional;
  CIT and the cohorts added in PR #11 are not included.
- The Erho directions come from a gene-level summary built outside the repository, after removing
  a technical intensity component. See modules/README.md.
