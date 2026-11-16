# Belfast / FASTMAN (GSE116918)

Pipeline for the Belfast prostate cancer cohort (Jain et al., *Ann Oncol* 2018, PMID 29045551):
diagnostic biopsies from men treated with radical radiotherapy and androgen-deprivation therapy (ADT),
profiled on the Almac Prostate Disease Specific Array. It is the only cohort in the project with a
**metastasis** endpoint and follow-up time, and the first to use one.

Code: `src/ds_bu_bioinformatics/belfast.py` (cohort adapter) and `belfast_pipeline.py` (builder).
Tests: `tests/test_belfast_pipeline.py` (synthetic data, no downloads).

## How to run

```bash
python -m ds_bu_bioinformatics.belfast_pipeline                  # snapshot belfast-v1 (downloads ~45 MB + ~400 MB)
python -m ds_bu_bioinformatics.belfast_pipeline --offline        # rebuild from cached, checksummed inputs
python -m ds_bu_bioinformatics.output_check --snapshot belfast-v1
```

The GPL25318 platform record is a single ~400 MB SOFT text file, so the first run takes a few minutes.
Inputs are cached with their URL, timestamp and SHA-256 in `data/raw/belfast-v1/manifest.json`;
offline rebuilds are byte-identical. Outputs use the shared layout (`data/processed/belfast-v1/`,
`data/anndata/belfast-v1/expression.h5ad`) and pass `check-outputs`.

## Sources

| Input | Source |
|---|---|
| Expression + clinical | GEO GSE116918 series matrix (248 samples) |
| Probe annotation | GEO GPL25318, Almac ADXPCv1a520642 (Ensembl release 75) |

Expression is the deposited RMA probeset-level log2 signal (range 0.65–14.42). No additional
normalization, transformation, imputation or averaging is applied.

## Samples

All 248 deposited profiles are primary prostate cancer diagnostic biopsies, one per patient
(GEO title `S1236-xxxxx`, used as `patient_id`). All 248 are kept, and every patient is eligible for
both endpoints. The pipeline fails if a profile has an unexpected tissue label or a duplicate patient.

## Probe annotation

The Almac array includes **antisense** probes, which measure transcripts from the opposite strand.
The annotation still lists a gene name for many of them (for example DYNLL1, PXDN), so mapping them
to that gene would mislabel antisense signal as gene expression. Only sense probes with exactly one
Entrez ID and one symbol are kept, and every probe keeps its own row.

| Status | Probes | Rule |
|---|---:|---|
| `mapped` | 43,152 | Sense (Fully Exonic or Includes Intronic), one Entrez ID, one symbol |
| `antisense_probe` | 7,836 | Orientation AntiSense, excluded even if a gene is named |
| `unannotated` | 5,696 | No gene symbol or Entrez ID |
| `multi_gene` | 4,776 | Symbol or Entrez ID contains `///` |
| `no_transcript_match` | 4,059 | No Genome Match / No Transcript Match |
| `gene_symbol_identity_conflict` | 12 | Same Entrez ID with different symbols, or the reverse |
| **Total deposited** | **65,531** | 43,152 rows kept, covering 16,174 unique genes |

Unknown orientation values stop the build. Every probe and its status are in `feature_mapping.csv.gz`.

## Endpoints

| Endpoint | Columns | Eligible | Events | Censored | Source fields |
|---|---|---:|---:|---:|---|
| **Metastasis (primary)** | `metastasis_*` | 248 | 22 | 226 | `met event (1=yes, 0=no)`, `follow-up time (met, months)` |
| BCR | `bcr_*` | 248 | 56 | 192 | `bcr event (1=yes, 0=no)`, `follow-up time (bcr, months)` |

Both follow the shared contract (`_event`, `_months`, `_eligible`, `_eligibility_reason`,
`_time_unit`). `metadata.json` records `primary_endpoint: metastasis`, so `module-scores` scores
Belfast against metastasis. BCR here is biochemical failure after radiotherapy + ADT, not
post-prostatectomy BCR.

The shared GEO builder (`geo_cohort_pipeline.py`) supports one BCR endpoint, so `belfast_pipeline.py`
reuses its feature selection, validation, CSV round trips and output protection unchanged, and adds the
metastasis endpoint and primary-endpoint metadata on top.

Clinical covariates are kept as deposited: age, PSA, clinical T stage (T1–T4; 25 Unknown) and total
Gleason score (6–10). `module-analysis` maps them to `age`, `gleason_8_plus` and `ct3_plus`
(clinical T3/T4).

## First results with the CampLab modules (metastasis endpoint)

130 of 166 modules scored (the rest have fewer than 5 measured genes). Adjusted models use age,
Gleason 8–10 and clinical T3+ on 223 complete cases with 22 events (5.5 events per parameter,
flagged as possibly unstable). The clinical-only model has a C-index of 0.61.

| Module | HR per SD (95% CI) | q | Adjusted HR (95% CI) | Adjusted q | C-index gain over clinical |
|---|---|---:|---|---:|---:|
| L58 · AR pathway | 0.39 (0.25–0.61) | 0.005 | 0.43 (0.27–0.69) | 0.068 | +0.15 |
| L60 · Neuroendocrine | 0.43 (0.28–0.66) | 0.006 | 0.48 (0.31–0.75) | 0.083 | +0.12 |
| L95 · TGFB pathway | 0.51 (0.33–0.79) | 0.11 | 0.57 (0.36–0.90) | 0.66 | +0.08 |
| L16 · Metallothionein | 0.59 (0.38–0.91) | 0.51 | 0.64 (0.41–1.02) | 0.98 | +0.06 |
| L12 · Cycling kinase | 1.14 (0.76–1.71) | — | 1.19 | — | — |
| L11 · DNA replication | 1.09 (0.70–1.68) | — | 0.98 | — | — |

- **AR pathway (L58) and L60 are strongly associated with lower metastasis risk** and keep most of
  their effect after adjustment (likelihood-ratio p = 0.0004 and 0.001 for adding the score to the
  clinical model). Both are FDR-significant unadjusted and close to it adjusted (q 0.07–0.08).
- **Proliferation modules (L12, L11) do not predict metastasis here**, unlike BCR/PFI in the other
  cohorts. This may reflect the small number of events, biopsy tissue, radiotherapy treatment, or a
  real difference between recurrence and metastasis biology.
- In the five-cohort meta-analysis (see [module_analysis.md](module_analysis.md)), L60 is protective
  in 5/5 cohorts (pooled HR 0.56) and L12 predicts higher risk in 5/5 (pooled HR 1.61; adjusted
  HR 1.35, 95% CI 1.13–1.59, adjusted q = 0.087).

## Caveats

- Only 22 metastasis events: adjusted estimates are unstable and wide.
- Diagnostic FFPE biopsy rather than prostatectomy tissue; treatment is radiotherapy + ADT.
- T stage is clinical, and Gleason is a total score without primary/secondary patterns, so the
  covariates are coarser than the pathological covariates in other cohorts.
- GEO event counts (BCR 56, metastasis 22) differ from Jain 2018 (65, 24); deposited values are used.
- The platform annotation is Ensembl release 75 (2013), so some symbols may not match current HGNC
  names used in the CampLab module list.
- No raw CEL quality control. C-indices are in-sample.