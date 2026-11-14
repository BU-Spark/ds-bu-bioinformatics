# TCGA-PRAD output column guide

Outputs of `uv run --locked tcga-prad` (snapshot `tcga-prad-v2`). The layout and naming follow the
Taylor guide in [data_dictionary.md](data_dictionary.md): the same file names, `feature_id`-keyed
expression, `patient__`/`sample__` prefixed source fields, and nullable endpoint columns built the
same way as Taylor's `dfs_*`. Field descriptions come from the cached cBioPortal attribute
definitions (`cbio_clinical_attributes.json`) and from the pipeline code.

Counts below come from the November 11, 2026 build.

## Reading the tables

- One row per selected sample in `samples.csv` and `clinical.csv`, in the same order as the expression columns.
- `patient__` = field from cBioPortal `prad_tcga_pan_can_atlas_2018` PATIENT records; `sample__` = SAMPLE records;
  `legacy_patient__` = Gleason fields from the TCGA Firehose Legacy study `prad_tcga`, matched on patient ID.
  Unprefixed fields are identifiers or are derived by this pipeline.
- Blank means missing/not supplied. It does not mean zero, negative, or event-free.
- Read identifiers as strings, for example `pd.read_csv(..., dtype={'sample_type_code': 'string', 'feature_id': 'string', 'gene_id': 'string'})`.
  Otherwise `01` becomes `1`.

## How TCGA differs from Taylor

| | Taylor (GSE21034) | TCGA-PRAD |
|---|---|---|
| Assay | Affymetrix exon array, RMA signal | RNA-seq, RSEM normalized counts (PanCanAtlas batch-corrected) |
| `feature_id` | Original GPL10264 probe ID | Entrez gene ID; `<entrez>_row<N>` when the ID repeats in the source file |
| Expression values | As deposited | log2(RSEM + 1) after clipping negative values to 0; raw RSEM in the h5ad `rsem` layer |
| Primary endpoint | DFS (`dfs_*`) | PFI (`pfi_*`); OS, DSS and DFI use the same column pattern |
| Audit of features | `feature_mapping.csv.gz` (RefSeq → Entrez mapping) | `feature_audit.csv` (no mapping step; records missing-value removals) |
| Sample selection | GEO `Primary tumor` label | Barcode type `01`, one per patient, then expression QC outliers removed |

## expression.csv.gz

19,088 feature rows × 484 sample columns. The index is `feature_id`; the columns are `sample_id`.

- Values are **log2(RSEM + 1)**, rounded to 4 decimals. The source values are PanCanAtlas batch-corrected
  RSEM counts. Correction leaves 40,717 small negative values (minimum −0.98); these are clipped to 0 before
  the log transform. The unclipped source values are kept in the AnnData layer `rsem`.
- 1,443 source rows have a missing value in at least one of 33 selected samples (batch-correction gaps).
  They are removed, not imputed, so every sample is kept. The filter only looks at the final selected samples,
  so a gap in an excluded sample does not remove a row. Sample QC (median expression, genes detected) uses
  the rows with no missing value across the primary tumors kept before QC.
- No rows are collapsed. Rows that share an Entrez ID or a symbol stay separate, as in Taylor:
  - 50 source rows share an Entrez ID with another row, for example 2258/FGF13 twice and 3192/HNRNPU vs NCRNA00201.
    Those rows hold different values (correlation 0.04–0.76), so all of them are kept. 47 of them remain
    after missing-value removal and use `<entrez>_row<N>` IDs.
  - 10 retained rows have no gene symbol.
- 19,064 unique Entrez genes.

## expression.h5ad (separate `data/anndata/<snapshot>/` folder)

| Component | Meaning |
|---|---|
| `X` | log2(RSEM + 1), samples × features, identical to `expression.csv.gz` |
| `layers['rsem']` | Source RSEM values for the same samples/features, unclipped |
| `obs` | `clinical.csv` fields indexed by `sample_id` |
| `var` | `genes.csv` fields indexed by `feature_id` |
| `uns['metadata_json']` | Source provenance, processing policies and schemas (same as `metadata.json`) |

## samples.csv

484 rows × 6 columns.

| Column | Meaning | Example |
|---|---|---|
| `sample_id` | TCGA sample barcode (15 characters); expression column key | TCGA-2A-A8VL-01 |
| `cbio_sample_id` | cBioPortal sample ID, checked against the sample inventory. Equal to `sample_id` in TCGA | TCGA-2A-A8VL-01 |
| `patient_id` | First 12 characters of the barcode; checked against the cBioPortal patient for that sample | TCGA-2A-A8VL |
| `sample_type_code` | Barcode characters 14–15. `01` = primary solid tumor | 01 |
| `platform_id` | Expression platform | IlluminaHiSeq_RNASeqV2 |
| `cbio_sample_type` | cBioPortal sample inventory type | Primary Solid Tumor |

## clinical.csv

484 rows × 89 columns. It starts with the identifiers `sample_id`, `cbio_sample_id`, `patient_id` and `sample_type_code`, as in `samples.csv`.

### Endpoints (derived)

Each endpoint has six columns: `<e>_event`, `<e>_months`, `<e>_eligibility_reason`, `<e>_eligible`,
`<e>_status` and `<e>_time_unit`. The decoding rules match Taylor DFS:

- Unknown status codes and negative or non-finite times stop the pipeline.
- A missing status or time stays missing; it is never treated as censored.
- A time of zero is kept but marked `zero_time_requires_review` and is not eligible.

| Prefix | Source fields | Event 1 / event 0 | Eligible | Events |
|---|---|---|---:|---:|
| `pfi` (**primary**) | `patient__PFS_STATUS` / `patient__PFS_MONTHS` | `1:PROGRESSION` / `0:CENSORED` | 484 | 91 |
| `os` | `patient__OS_STATUS` / `patient__OS_MONTHS` | `1:DECEASED` / `0:LIVING` | 484 | 10 |
| `dss` | `patient__DSS_STATUS` / `patient__DSS_MONTHS` | `1:DEAD WITH TUMOR` / `0:ALIVE OR DEAD TUMOR FREE` | 482 | 5 |
| `dfi` | `patient__DFS_STATUS` / `patient__DFS_MONTHS` | `1:Recurred/Progressed` / `0:DiseaseFree` | 330 | 30 |

| Column | Meaning |
|---|---|
| `<e>_event` | 1 = event, 0 = censored at `<e>_months`; blank = status missing |
| `<e>_months` | Numeric time from the source field; no unit conversion |
| `<e>_eligibility_reason` | `eligible`, `missing_status`, `missing_time` or `zero_time_requires_review` |
| `<e>_eligible` | True only when the reason is `eligible` |
| `<e>_status` | Source status text without the numeric prefix |
| `<e>_time_unit` | `months` |

These are the TCGA Pan-Cancer Clinical Data Resource endpoints as cBioPortal exposes them. cBioPortal's
`PFS_*` fields are the CDR progression-free interval (PFI), and `DFS_*` is the disease-free interval (DFI),
which is incomplete in PRAD. PFI is the client-agreed endpoint. Do not relabel it BCR or DFS when comparing
with other cohorts.

### Covariates (derived)

| Column | Meaning | Missing |
|---|---|---:|
| `age` | Numeric `patient__AGE` (years at diagnosis). Missing tokens become blank; any other non-numeric value or one outside 0–120 stops the build | 0 |
| `path_t_group` | T category (T2/T3/T4) taken from `patient__PATH_T_STAGE` | 7 |
| `gleason_score` | Integer `legacy_patient__GLEASON_SCORE`. Missing tokens become blank; any other non-integer value or one outside 2–10 stops the build | 0 |

### Gleason (TCGA Firehose Legacy, `prad_tcga`)

| Column | Meaning | Missing |
|---|---|---:|
| `legacy_patient__GLEASON_SCORE` | Total Gleason score | 0 |
| `legacy_patient__GLEASON_PATTERN_PRIMARY` | Primary Gleason pattern | 0 |
| `legacy_patient__GLEASON_PATTERN_SECONDARY` | Secondary Gleason pattern | 0 |
| `legacy_patient__GLEASON_PATTERN_TERTIARY` | Tertiary pattern; blank usually means no tertiary pattern was reported | 335 |

### Source fields

Values are kept verbatim. The descriptions are cBioPortal's.

| Column | Description | Missing |
|---|---|---:|
| `patient__AGE` | Age at which the condition was first diagnosed | 0 |
| `patient__PATH_T_STAGE` | Pathologic T stage (AJCC) | 7 |
| `patient__PATH_N_STAGE` | Pathologic N stage (AJCC) | 71 |
| `patient__PFS_STATUS`, `patient__PFS_MONTHS` | Progression-free status / months | 0 |
| `patient__OS_STATUS`, `patient__OS_MONTHS` | Overall survival status / months since initial diagnosis | 0 |
| `patient__DSS_STATUS`, `patient__DSS_MONTHS` | Disease-specific survival status / months | 2 / 0 |
| `patient__DFS_STATUS`, `patient__DFS_MONTHS` | Disease-free status / months since initial treatment | 154 |
| `patient__NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT` | New tumor event after initial treatment (Yes/No) | 57 |
| `patient__PERSON_NEOPLASM_CANCER_STATUS` | Tumor Free / With Tumor at last contact | 85 |
| `patient__RADIATION_THERAPY`, `patient__HISTORY_NEOADJUVANT_TRTYN` | Treatment history | 45 / 0 |
| `patient__BUFFA_HYPOXIA_SCORE`, `patient__RAGNUM_HYPOXIA_SCORE`, `patient__WINTER_HYPOXIA_SCORE` | mRNA-based hypoxia scores | 154 |
| `patient__SUBTYPE` | TCGA molecular subtype label | 13 |
| `patient__RACE`, `patient__ETHNICITY`, `patient__GENETIC_ANCESTRY_LABEL` | Demographics / genetic ancestry | 328 / 332 / 15 |
| `patient__DAYS_LAST_FOLLOWUP`, `patient__DAYS_TO_BIRTH`, `patient__DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS` | Day offsets relative to diagnosis | 8 / 11 / 31 |
| other `patient__*` | ICD codes, consent, form date, PanCan freeze flag, sample count, sex, other IDs | 0 |
| `sample__SAMPLE_TYPE` | cBioPortal sample type; must be `Primary` for selected samples | 0 |
| `sample__TUMOR_TYPE` | Histologic tumor type (acinar, other subtypes) | 0 |
| `sample__FRACTION_GENOME_ALTERED`, `sample__ANEUPLOIDY_SCORE`, `sample__TBL_SCORE` | Copy-number burden measures | 4 / 21 / 2 |
| `sample__MUTATION_COUNT`, `sample__TMB_NONSYNONYMOUS` | Mutation burden | 0 |
| `sample__MSI_SCORE_MANTIS`, `sample__MSI_SENSOR_SCORE` | Microsatellite instability scores | 0 |
| other `sample__*` | Cancer type labels, OncoTree code, somatic status, tissue source site, collection indicators | 0–2 |

## genes.csv

One row per expression row, in the same order. 19,088 rows × 6 columns.

| Column | Meaning |
|---|---|
| `feature_id` | Expression row key (see above) |
| `gene_id` | Entrez gene ID (string); can repeat |
| `gene_symbol` | Source `Hugo_Symbol`; can repeat or be blank. Not re-validated against HGNC |
| `source_row` | 1-based data row in `data_mrna_seq_v2_rsem.txt` |
| `mean_log2` | Mean log2(RSEM + 1) across the selected samples |
| `pct_samples_expressed` | % of selected samples with log2(RSEM + 1) > 1 |

## feature_audit.csv

Every source row (20,531): `feature_id`, `source_row`, `gene_id`, `gene_symbol`, `missing_samples`
(the number of samples with a missing value), `selected`, and `reason` (`selected` or `missing_values`).

## sample_audit.csv

Every expression column (493): `sample_id`, `patient_id`, `sample_type_code`, `selected`, `selection_reason`,
`median_log2`, `genes_detected`, `cbio_sample_id`, `cbio_sample_type`.

`selection_reason` takes these values:

- `primary_tumor` (selected; 484)
- `qc_outlier` (9)
- `not_primary_tumor`
- `extra_primary_for_patient`
- `not_in_cbio_inventory`
- `no_clinical_record`

The last four do not occur in this build. A QC outlier lies more than 3 SD from the mean of the
candidate samples in median log2 expression or in genes detected (log2(RSEM + 1) > 1).

## module_coverage.csv (only with `--module-file`)

One row per module: `genes`, `found` (symbol present in `genes.csv`), `lost_to_missing_values`
(symbol only in removed rows), and `pct_found`.

## validation.json, metadata.json, qc_summary.md, checksums.sha256

- `validation.json` holds the structural checks (alignment, unique IDs, finite values, valid event codes), plus:
  - counts per endpoint (eligible records only; `counts.pfi_events`/`pfi_censored` equal
    `endpoints.pfi.events`/`censored`)
  - feature and sample dispositions
  - clinical missingness
  - export and AnnData round-trip results
- `metadata.json` holds:
  - the source URLs (the RSEM file is pinned to cBioPortal datahub commit `86690e1`) and the raw-manifest hash
  - the expression transform
  - the endpoint code maps and cBioPortal definitions
  - sample policy, field provenance and limitations
  - output schemas and code hashes
- `qc_summary.md` is a readable summary of the same counts.
- `checksums.sha256` holds SHA-256 hashes of every file in the folder.
