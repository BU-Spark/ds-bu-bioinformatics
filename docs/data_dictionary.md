# Taylor/MSKCC output column guide

This guide describes the actual `taylor-v2` output headers and pipeline behavior. Source descriptions come from the locally cached cBioPortal attribute definitions in `metadata.json`, GEO fields, and the implementation in `src/ds_bu_bioinformatics/`. It does not invent definitions absent from those sources.

## Reading the tables

- A row is one item: a mapped feature in expression/genes, a selected profile in clinical/samples, a source profile in the sample audit, or a feature/candidate record in the mapping audit.
- `patient__` means a field downloaded from cBioPortal patient records. `sample__` means a field downloaded from cBioPortal sample records. Unprefixed source fields are generally GEO fields; `dfs_*` and selection fields are derived by this pipeline.
- Blank means missing/not supplied, not zero, normal, negative, or disease-free. Source definitions may be incomplete.
- QC means quality control. Here it covers downloaded file integrity, identity matching, valid values, alignment and export checks. It does not include raw CEL assay-quality testing.
- DFS means disease-free survival. Source-defined observation time is months since initial treatment. Exact recurrence criteria, death handling and a surgery-specific origin are not fully verified. Event 0 means no event recorded through the observation time, not permanent cure.

## expression.csv.gz

The canonical expression CSV has one row per retained GPL10264 feature, indexed
by `feature_id`, and one measurement column per selected sample (131 in the frozen
`taylor-v2` cohort). Its shape is recorded in `validation.json` as
`retained_feature_rows` × `counts.primary_profiles`. Join the feature ID to
`genes.csv.feature_id`; gene IDs and symbols may repeat across rows. Each value
is the original deposited signal for that feature/sample; no averaging is applied.

GEO reports RMA background adjustment and quantile normalization. The exact
deposited log base is unverified. No log transform or additional normalization is
applied. Values are not read counts or percentages.

## expression.h5ad (separate data/anndata/<snapshot>/ folder)

| Component | Meaning |
|---|---|
| `X` | Expression transposed to samples × retained features, preserving numeric values. |
| `obs` | Clinical/sample fields from clinical.csv, indexed by sample_id in expression order. Missing clinical strings remain missing; nullable floating times use NaN. |
| `var` | Gene ID, symbol and original/query RefSeq annotations, indexed by unique feature_id. Gene IDs and symbols may repeat. |
| `uns['metadata_json']` | JSON string containing source provenance, schemas and processing policies. |

## samples.csv

131 data rows × 6 columns.

| Column | Meaning | Example | Missing rows |
|---|---|---|---:|
| `sample_id` | Study sample identifier and the key used for expression columns. Example: PCA0001. It is a sample key even where its text equals the patient ID. | PCA0001; PCA0002 | 0 |
| `geo_accession` | GEO sample accession (GSM), linking the profile to its GEO source record. | GSM527858; GSM527859 | 0 |
| `cbio_sample_id` | Actual cBioPortal sample identifier, matched exactly to the study sample ID and checked against the sample inventory. | PCA0001; PCA0002 | 0 |
| `patient_id` | Actual cBioPortal patient identifier. Distinct concept from sample_id. The selected dataset has one sample per patient. | PCA0001; PCA0002 | 0 |
| `platform_id` | Expression platform accession. GPL10264 in this dataset. | GPL10264 | 0 |
| `tumor_type` | Original GEO tumor-type label. Only Primary tumor profiles enter the final selected dataset. | Primary tumor | 0 |

## clinical.csv

131 data rows × 46 columns.

| Column | Meaning | Example | Missing rows |
|---|---|---|---:|
| `sample_id` | Study sample identifier and the key used for expression columns. Example: PCA0001. It is a sample key even where its text equals the patient ID. | PCA0001; PCA0002 | 0 |
| `geo_accession` | GEO sample accession (GSM), linking the profile to its GEO source record. | GSM527858; GSM527859 | 0 |
| `title` | Original GEO sample title, retained as descriptive source text. | Prostate tumor PCA0001 transcript; Prostate tumor PCA0002 transcript | 0 |
| `platform_id` | Expression platform accession. GPL10264 in this dataset. | GPL10264 | 0 |
| `tissue` | Original GEO tissue label. | Prostate | 0 |
| `disease_status` | Original GEO disease-status description. This is not the derived DFS recurrence event. | prostate cancer | 0 |
| `tumor_type` | Original GEO tumor-type label. Only Primary tumor profiles enter the final selected dataset. | Primary tumor | 0 |
| `biopsy_gleason_grade` | GEO-reported biopsy Gleason grade/score field. Retained separately from the cBioPortal surgical-specimen score. Do not assume the two measurements are interchangeable. | 6; 7 | 0 |
| `clint_stage` | Original GEO clinical-stage field, with its original spelling retained. Do not automatically replace the separate cBioPortal clinical-stage field with this value. | T2A; T1C | 0 |
| `pathological_stage` | Original GEO pathological-stage field. Retained separately from patient__PATH_T_STAGE. | T3C; T2A | 0 |
| `cell_line` | Original GEO cell-line metadata field. Entirely missing among the selected primary tumors. Blank does not establish a measured biological value. |  | 131 |
| `cell_type` | Original GEO cell-type metadata field. Entirely missing among selected primary tumors. |  | 131 |
| `selected` | Pipeline selection flag: True when GEO labels the profile Primary tumor. Final clinical rows are selected; the audit also includes excluded rows. | True | 0 |
| `selection_reason` | Pipeline sample-selection explanation: primary_tumor, not_primary_tumor, or missing_tumor_type. | primary_tumor | 0 |
| `cbio_sample_id` | Actual cBioPortal sample identifier, matched exactly to the study sample ID and checked against the sample inventory. | PCA0001; PCA0002 | 0 |
| `sample__CANCER_TYPE` | cBioPortal broad cancer-type label. | Prostate Cancer | 0 |
| `sample__CANCER_TYPE_DETAILED` | cBioPortal more specific cancer diagnosis label. | Prostate Adenocarcinoma | 0 |
| `sample__COMPLETE_DATA` | Original cBioPortal Complete Data flag. Source metadata gives no detailed completeness rule. This is not the pipeline QC result. | Yes; No | 1 |
| `sample__ERG_FUSION_ACGH` | Source ERG-fusion call labeled ACGH. Preserve the source Positive/Negative values. Detailed calling criteria are not defined in the cached attribute description. | Negative; Positive | 0 |
| `sample__ERG_FUSION_GEX` | Source ERG-fusion call labeled GEX (gene expression). Preserve separately from the ACGH call. Detailed calling criteria are not defined in the cached attribute description. | Negative; Positive | 0 |
| `sample__FRACTION_GENOME_ALTERED` | Source fraction of the genome reported as altered. A fraction such as 0.18 represents about 18%; the cached definition does not establish the detailed measurement/calling method. | 0.1828; 0.0761 | 22 |
| `sample__GLEASON_SCORE` | cBioPortal radical-prostatectomy Gleason score, as labeled by the source metadata. Separate from the GEO biopsy field. | 7; 8 | 1 |
| `sample__GLEASON_SCORE_1` | Source additional Gleason score field ending _1. Values are retained verbatim. Its cached definition repeats the overall Gleason definition and does not explicitly verify that it is the primary pattern. | 4; 5 | 2 |
| `sample__GLEASON_SCORE_2` | Source additional Gleason score field ending _2. Values are retained verbatim. Its cached definition does not explicitly verify that it is the secondary pattern. | 3; 4 | 1 |
| `sample__MUTATION_COUNT` | Source-reported mutation count for the sample. Blank is unknown/unavailable, not zero. The exact mutation filtering scope is not specified by the cached definition. | 6; 1 | 101 |
| `sample__ONCOTREE_CODE` | Source OncoTree cancer-classification code. PRAD is the observed code in this selected dataset. | PRAD | 0 |
| `sample__SAMPLE_CLASS` | Source sample class, such as Tumor, cell line, or xenograft. Distinct from primary versus metastatic status. | Tumor | 0 |
| `sample__SAMPLE_TYPE` | Source sample type, such as Primary, metastatic, or normal. Selected records must be Primary when this field is supplied. | Primary | 0 |
| `sample__SEQUENCED` | Source flag indicating sequencing status. This does not imply the expression assay itself is RNA sequencing; this pipeline uses a processed microarray matrix. | Yes; No | 0 |
| `sample__SEQUENCING` | Source text describing sequencing methods. Preserved as supplied, without inferring whole-genome or whole-exome coverage. | Sanger and Sequenom; Sequenom | 0 |
| `sample__SOMATIC_STATUS` | Original source somatic-status label, with values such as Matched or Unmatched. Detailed matching criteria are not provided in the cached attribute description. | Unmatched; Matched | 0 |
| `sample__TMB_NONSYNONYMOUS` | Source nonsynonymous tumor mutational burden value. Exact units and filtering method are not specified by the cached attribute definition; do not assume or relabel units. | 0; 15.74865152 | 56 |
| `patient_id` | Actual cBioPortal patient identifier. Distinct concept from sample_id. The selected dataset has one sample per patient. | PCA0001; PCA0002 | 0 |
| `patient__CLIN_T_STAGE` | Patient-level cBioPortal clinical primary-tumor T stage assessed before treatment, according to the source attribute definition. | T2A; T1C | 0 |
| `patient__COPY_NUMBER_CLUSTER` | Source copy-number cluster assignment. The cluster number is a group label, not a quantitative severity score. | 6; 2 | 0 |
| `patient__DFS_MONTHS` | Original source DFS time in months since initial treatment. Exact treatment/surgery origin and detailed recurrence rules remain unverified. | 18.5; 58.02 | 0 |
| `patient__DFS_STATUS` | Original source DFS status string. 0:DiseaseFree maps to event 0; 1:Recurred or 1:Recurred/Progressed maps to event 1. | 1:Recurred; 0:DiseaseFree | 0 |
| `patient__PATH_T_STAGE` | Patient-level cBioPortal pathological primary-tumor T-stage code. Kept separately from the original GEO pathological_stage. | T3C; T2A | 0 |
| `patient__SAMPLE_COUNT` | Source number of samples per patient in the cBioPortal study. Distinct from the pipeline-selected sample count. | 1 | 0 |
| `patient__SEX` | Source-recorded sex. | Male | 0 |
| `dfs_event` | Derived numeric DFS event: 1 = recurrence/progression recorded by source; 0 = source DiseaseFree, censored at recorded time. Missing status remains missing and is never converted to 0. | 1; 0 | 0 |
| `dfs_months` | Numeric version of patient__DFS_MONTHS. No unit conversion. For an event, the recorded time to event; for event 0, the recorded time to censoring. | 18.5; 58.020000000000003 | 0 |
| `dfs_eligibility_reason` | Pipeline explanation of whether DFS fields can be used: eligible, missing_status, missing_time, or zero_time_requires_review. | eligible | 0 |
| `dfs_eligible` | True only when dfs_eligibility_reason is eligible. This means complete usable status/time under the pipeline rules, not eligibility for every possible analysis. | True | 0 |
| `dfs_status` | Source status text after removing the numeric prefix, for example Recurred or DiseaseFree. Derived from patient__DFS_STATUS. | Recurred; DiseaseFree | 0 |
| `dfs_time_unit` | Explicit unit for dfs_months: months. | months | 0 |

## genes.csv

One row per retained expression feature, in the same order; five columns.

| Column | Meaning |
|---|---|
| `feature_id` | Unique original GPL10264 feature ID; expression row key. |
| `gene_id` | Entrez gene ID annotation, stored as a string; may repeat. |
| `gene_symbol` | MyGene gene symbol annotation; may repeat. Not claimed HGNC-approved. |
| `refseq_original` | Original GPL10264 RefSeq accession. |
| `refseq_query` | Lookup accession after removing only a terminal numeric version. |

## sample_audit.csv

185 data rows × 16 columns.

| Column | Meaning | Example | Missing rows |
|---|---|---|---:|
| `geo_accession` | GEO sample accession (GSM), linking the profile to its GEO source record. | GSM527858; GSM527859 | 0 |
| `title` | Original GEO sample title, retained as descriptive source text. | Prostate tumor PCA0001 transcript; Prostate tumor PCA0002 transcript | 0 |
| `platform_id` | Expression platform accession. GPL10264 in this dataset. | GPL10264 | 0 |
| `sample_id` | Study sample identifier and the key used for expression columns. Example: PCA0001. It is a sample key even where its text equals the patient ID. | PCA0001; PCA0002 | 0 |
| `tissue` | Original GEO tissue label. | Prostate | 0 |
| `disease_status` | Original GEO disease-status description. This is not the derived DFS recurrence event. | prostate cancer; normal adjacent benign prostate | 6 |
| `tumor_type` | Original GEO tumor-type label. Only Primary tumor profiles enter the final selected dataset. | Primary tumor; Metastasis | 35 |
| `biopsy_gleason_grade` | GEO-reported biopsy Gleason grade/score field. Retained separately from the cBioPortal surgical-specimen score. Do not assume the two measurements are interchangeable. | 6; 7 | 36 |
| `clint_stage` | Original GEO clinical-stage field, with its original spelling retained. Do not automatically replace the separate cBioPortal clinical-stage field with this value. | T2A; T1C | 40 |
| `pathological_stage` | Original GEO pathological-stage field. Retained separately from patient__PATH_T_STAGE. | T3C; T2A | 44 |
| `cell_line` | Original GEO cell-line metadata field. Entirely missing among the selected primary tumors. Blank does not establish a measured biological value. | DU145; LNCaP | 181 |
| `cell_type` | Original GEO cell-type metadata field. Entirely missing among selected primary tumors. | Xenograft LNCaP104R; Xenograft LNCaP104S | 183 |
| `selected` | Pipeline selection flag: True when GEO labels the profile Primary tumor. Final clinical rows are selected; the audit also includes excluded rows. | True; False | 0 |
| `selection_reason` | Pipeline sample-selection explanation: primary_tumor, not_primary_tumor, or missing_tumor_type. | primary_tumor; not_primary_tumor | 0 |
| `cbio_sample_id` | Actual cBioPortal sample identifier, matched exactly to the study sample ID and checked against the sample inventory. | PCA0001; PCA0002 | 29 |
| `patient_id` | Actual cBioPortal patient identifier. Distinct concept from sample_id. The selected dataset has one sample per patient. | PCA0001; PCA0002 | 29 |

## feature_mapping.csv.gz

43,419 data rows × 12 columns.

| Column | Meaning | Example | Missing rows |
|---|---|---|---:|
| `feature_id` | Original numeric GPL10264 feature identifier before mapping selection. It is not an Entrez gene ID. | 1; 10 | 0 |
| `refseq_original` | Original RefSeq accession supplied in the GPL10264 GB_ACC annotation field. | NM_000014; NM_000023 | 0 |
| `refseq_query` | RefSeq accession used for lookup after removing only a final numeric version suffix, if present. | NM_000014; NM_000023 | 0 |
| `candidate_gene_id` | Entrez ID proposed by a returned MyGene candidate. It becomes a final mapping only when selected is True. | 2; 6442 | 17,519 |
| `candidate_symbol` | Gene symbol proposed by the candidate. | A2M; SGCA | 17,519 |
| `candidate_taxid` | Candidate species taxonomy identifier. Human candidates must be 9606. | 9606 | 17,519 |
| `candidate_score` | MyGene search relevance score. Not expression, probability, or a confidence threshold. The pipeline does not simply select the highest score. | 10.910185999999999; 19.682576999999998 | 17,519 |
| `candidate_refseq_rna` | Semicolon-separated RefSeq RNA accessions returned for this candidate, used to verify exact version-stripped accession agreement. | NM_000014.6;NM_001347423.2;NM_001347424.2;NM_001347425.2;XM_006719056.4;XM_054371391.1; NM_000023.4;NM_001135697.3;NR... | 17,519 |
| `mapping_status` | Mapping disposition, such as mapped, unmapped, invalid_refseq, ambiguous, invalid_candidate, refseq_mismatch, or gene_symbol_identity_conflict. | mapped; unmapped | 0 |
| `selected` | True only when this feature/candidate mapping passes the explicit identity rules and is retained in the expression matrix. Different from sample selection. | True; False | 0 |
| `reason` | Explanation of the mapping decision; currently mirrors mapping_status. | mapped; unmapped | 0 |
| `mapping_source` | Cached service response used for the mapping: MyGene.info v3. | MyGene.info v3 cached response | 0 |

## validation.json: the QC summary

This is a JSON report rather than a CSV table. `clinical_missingness` lists a missing-value count for every clinical column.

| Field | Meaning |
|---|---|
| `structurally_valid` | Alignment, unique identifiers, finite expression, selection and DFS validity checks passed. |
| `counts` | Final selected profiles, unique Entrez genes (gene_rows) and source-coded event/non-event counts. |
| `retained_feature_rows` | Number of separate mapped feature rows in CSV / variables in AnnData. |
| `anndata_round_trip_valid` | H5AD reload preserved values, sample/feature order, clinical annotations and provenance. |
| `geo_profiles` | All source GEO expression profiles before selection: 185. |
| `source_features` | Original measured platform features before mapping selection: 43,419. |
| `feature_dispositions` | Counts of distinct features by mapping status, not candidate-row counts. |
| `eligible_dfs_profiles` | Number with DFS status and positive recorded DFS time under current rules. |
| `clinical_missingness` | Count of blank/missing records per clinical column. |
| `expression_range` | Minimum and maximum final expression measurement. |
| `identical_patient_records_collapsed` / `identical_sample_records_collapsed` | Identical repeated source clinical records removed. Conflicting duplicates fail. |
| `export_round_trip_valid` | Reloading the exported files preserved the checked data. |
| `legacy_baseline` | Earlier output counts, retained for comparison only. |
| `baseline_differences` | Current minus earlier counts. No target count is forced. |
| `expression_scale_verified` | False: exact deposited expression scale/log base remains unresolved. |
| `precise_dfs_definition_verified` | False: precise recurrence criteria are not fully verified. |

## Supporting provenance files

- `metadata.json` stores the source attribute definitions, field provenance, output schemas, processing policies, limitations, code hashes and raw-manifest hash.
- `checksums.sha256` contains the expected SHA-256 digest and filename of each tracked final output. It detects byte changes. It does not establish biological validity. The AnnData folder has its own checksum manifest.
- `data/raw/taylor-v2/manifest.json` is the **download manifest**. It records source requests, URLs, retrieval timestamps and input hashes. It is distinct from the **sample manifest**, `samples.csv`.

## What to point to during the demo

1. In `samples.csv`, follow PCA0001 to GSM527858 and patient PCA0001.
2. In `clinical.csv`, show `patient__DFS_STATUS = 1:Recurred`, `dfs_event = 1`, and `dfs_months = 18.5` for that profile.
3. In `expression.csv.gz`, show the PCA0001 column. Its values describe gene expression, not the clinical outcome.
4. In `genes.csv`, connect each expression `feature_id` to its gene ID, symbol and RefSeq annotation. Show repeated symbols on separate rows.
5. In `validation.json`, show the counts, missingness and alignment result. State the remaining expression-scale and DFS-definition limits.
