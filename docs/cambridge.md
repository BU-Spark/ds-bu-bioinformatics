# Cambridge / CamCap GSE70768

This adapter downloads the processed GSE70768 Series Matrix, GPL10558
Illumina HumanHT-12 v4 annotation, and the original publication clinical supplement
(PMC4588396, Supplementary Table 2). It selects primary tumors independently of
outcome completeness, preserves separate probes sharing a gene annotation, and
exports aligned CSV and H5AD datasets. It does not perform raw assay processing.

## Run

From a checkout containing the Cambridge adapter:

```bash
uv sync --locked
uv run --locked cambridge-gse70768
```

Default snapshot: `cambridge-v2`. The output structure parallels Taylor:

```text
data/raw/cambridge-v2/         Series Matrix, platform SOFT, publication ZIP, manifest.json
data/interim/cambridge-v2/     All source features for the selected primary tumors
data/processed/cambridge-v2/   expression.csv.gz, genes.csv, clinical.csv,
                              samples.csv, feature_mapping.csv.gz,
                              sample_audit.csv, metadata.json, validation.json,
                              checksums.sha256
data/anndata/cambridge-v2/     expression.h5ad, metadata.json, checksums.sha256
```

As for Taylor, H5AD uses samples × features in `X`, clinical/sample fields in
`obs`, probe/gene annotations in `var`, and JSON provenance in
`uns['metadata_json']`. Original ILMN probe IDs are unique variable/CSV row keys;
Entrez IDs and symbols may repeat. Deposited values are not averaged or imputed.
The pipeline validates reloaded CSV and H5AD exports before publishing any format.
Existing differing outputs are protected.

After acquisition, reproduce independently without network access:

```bash
uv run --offline --locked cambridge-gse70768 --offline \
  --output-root /tmp/cambridge-check/processed \
  --interim-root /tmp/cambridge-check/interim \
  --anndata-root /tmp/cambridge-check/anndata
```

Copy the entire `data/raw/cambridge-v2/` cache and manifest to another machine for
an exact-input rebuild. The committed/reference manifest alone does not contain
source payloads. Use `--snapshot cambridge-v3` for a deliberate source refresh.
A fresh online acquisition may retrieve changed upstream annotation. Use new
output roots when code/metadata changes rather than replacing frozen exports.

## Selection and clinical fields

- Select GEO sample type `Tumour` **and** title prefix
  `tumour tissue_robotic radical prostatectomy_`.
- Exclude the CRPC title prefix and benign specimens. Known benign titles use
  either `matched benign tissue_...` or `benign tissue_...`; other title/type
  combinations fail for explicit review.
- Keep every selected primary, even when BCR status/time is missing. The
  `bcr_eligible` flag defines a downstream analysis subset; it does not select
  expression columns.
- `sample_id` equals the GEO accession. `patient_id` is the cohort subject label
  parsed from the final underscore-separated title segment and matched to the
  publication SubjectID prefix. It is not an external registry patient ID. Reject
  repeated primary GEO subject IDs or conflicting publication clinical records.
  Collapse multiple publication tumor-region rows only when their clinical values
  agree exactly, and preserve every contributing publication row ID.
- Preserve original characteristics with `geo__` prefixes. Characteristic keys
  are parsed per sample because GEO characteristic rows can contain different
  keys in different columns.
- Preserve `geo__` characteristics unchanged, but use the publication's BCR
  status, event time and follow-up **together for every primary tumor**. The
  publication has later follow-up and higher event-time precision than GEO.
- `paper__subject_id` preserves matched original publication row IDs;
  `paper__bcr_status`, `paper__bcr_type`, `paper__bcr_months` and
  `paper__followup_months` preserve the selected clinical source fields.
- Map publication BCR `Y → 1`, `N → 0`; unknown nonmissing codes fail. Event
  duration comes from `paper__bcr_months`; censoring comes from
  `paper__followup_months`. Both are source months. `bcr_clinical_source` identifies
  the source for each observation. No event times are rounded.
- Missing status/time remains nullable. Zero time is retained but ineligible.
  Negative, nonfinite or malformed times fail. An event later than recorded
  follow-up is retained with reason
  `event_after_recorded_followup_requires_review` and is ineligible.

The source-coded endpoint is **BCR**, distinct from Taylor's DFS. Publication BCR
subtypes include early relapse and treatment-triggered relapse. The precise
uniform PSA threshold, confirmation rule and death handling are not established
by this adapter; use the preserved BCR subtype field for source review. The
clinical source is the paper's complete table, rather than a mixture of selected
GEO and paper values. The workbook is read from the checksummed source ZIP using
its cached cell values; no spreadsheet recalculation or editing is performed.

## Mapping and missing expression

The frozen GPL10558 columns `Species`, `Entrez_Gene_ID`, `Symbol` and `RefSeq_ID`
provide annotation. Require human species, a single numeric Entrez ID and a single
nonempty symbol. Quarantine gene IDs with inconsistent symbols, or symbols shared
across different gene IDs. Features sharing the same gene ID/symbol remain separate.
RefSeq is preserved as deposited, not used to assert a newly validated identity.

This differs from Taylor's frozen MyGene/exact-RefSeq mapping. Cambridge does not
refresh annotations through MyGene, and deposited labels are not claimed current
or HGNC-approved. Mapping dispositions and the source hash are recorded, so a
future annotation refresh can be compared explicitly.

Reject infinities and malformed numeric cells. Recognized missing expression is
allowed only in source/intermediate data. Exclude an entire feature if it has any
missing value across **all selected primary tumors**, including those with missing
clinical outcomes. Do not exclude samples, fill cells, average probes, normalize,
transform or adjust batches. Source text reports log2 transformation, quantile
normalization and ComBat; the adapter adds none of these operations.

`feature_mapping.csv.gz` audits every deposited matrix feature, including:

| Field | Meaning |
|---|---|
| `feature_id` | Original ILMN probe ID, in source order |
| `gene_id`, `gene_symbol`, `refseq_original`, `species` | Deposited GPL10558 annotation |
| `mapping_status` | mapped, invalid_or_ambiguous_annotation, or gene_symbol_identity_conflict |
| `missing_selected_values` | Missing cells across all selected primary tumors |
| `expression_status` | complete or missing_selected_expression |
| `selected` | Passes annotation policy and has complete selected expression |
| `mapping_source` | Explicit deposited-platform source; no MyGene refresh |

## Observed source issues and counts

The November 6, 2026 acquisition contains 199 profiles: **112 primary tumors,
13 CRPC and 74 benign**. Primary profiles have 19 source-coded BCR events and
93 non-events. One non-event has no follow-up. Thus **111 have complete positive
time**, but four events have a time later than their recorded total follow-up:

| GEO accession | Subject | BCR months | Recorded follow-up months |
|---|---|---:|---:|
| GSM1817747 | TB08.1026 | 57 | 38.8931506849315 |
| GSM1817773 | TB09.1408 | 31 | 30.3780821917808 |
| GSM1817810 | TB11.2182 | 21 | 14.6301369863014 |
| GSM1817820 | TB12.0560 | 19 | 2.03835616438356 |

The publication resolves all four timing conflicts and supplies follow-up for
GSM1817832 / TB12.2336, which lacks it in GEO. Its corresponding tumor rows are
`TB08.1026_17_T`, `TB09.1408_10_T`, `TB11.2182_1_T`, `TB12.0560_9_T` and
`TB12.2336_12_T`. The paper follow-up values are 57.369863013698634,
45.5013698630137, 20.71232876712329, 19.463013698630135 and
10.454794520547946 months, respectively. The paper also reports event times at
higher precision. All original GEO values remain visible for comparison.

With the complete publication clinical source, **all 112 primary tumors are
BCR-eligible: 19 events and 93 censored**. All 112 are retained in expression,
clinical CSV and H5AD. The previously reported 111/19 figures described GEO-only
completeness, not the subsequently reconciled publication table. Counts are
observations, not selection targets.

Snapshot `cambridge-v1` was the initial GEO-only diagnostic build. The default
`cambridge-v2` snapshot adds the publication source and supersedes that clinical
build. Its generated outputs were moved to the ignored
`.codex-build/cambridge-geo-only/` archive; active Cambridge outputs use v2.
Taylor snapshots remain separate.

There are 47,323 deposited probes. Across the 112 primary tumors, 12 probes have
48 missing cells; 11 otherwise-mapped probes are excluded for missing expression.
The annotation rules retain **41,610 separate features representing 29,604 unique
Entrez IDs** after missing-expression exclusion. These counts differ from the
exploratory check because production additionally enforces human species and
consistent gene/symbol identity. Counts are observations, not validation targets.

## Verification

```bash
uv run --locked python -m unittest discover -s tests -v
```

Synthetic Cambridge tests cover independent specimen selection, missing clinical
outcomes, event/censoring rules, conflicting clinical times, missing probes in
ineligible primaries, missing values only in excluded specimens, ambiguous mapping,
separate same-gene probes, malformed sources, CSV/H5AD reloads, offline reproduction,
cohort cache identity, publication subject joins/conflicting clinical records,
and protection of all output formats. Existing Taylor tests
also run. Snapshot counts are reported rather than hardcoded in tests.

Sources: [GSE70768](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70768),
[GPL10558](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558),
[primary publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC4588396/).
