# Prostate cancer expression and clinical pipelines

Reproducible acquisition, annotation and clinical matching for prostate cancer
cohorts. Each pipeline selects primary tumors, checks sample-to-patient identities,
preserves separate features sharing a gene, and exports aligned CSV and AnnData
(`.h5ad`) datasets with the same file layout.

| Cohort | Command | Sources | Guide |
|---|---|---|---|
| Taylor/MSKCC | `taylor-gse21034` | GEO **GSE21034 / GPL10264**, cBioPortal **prad_mskcc** | [data dictionary](docs/data_dictionary.md) |
| Sboner/Swedish WW | `sboner-gse16560` | GEO **GSE16560 / GPL5474** (clinical from GEO sample characteristics) | [Sboner / GSE16560](#sboner--gse16560-gpl5474) |
| TCGA-PRAD | `tcga-prad` | cBioPortal **prad_tcga_pan_can_atlas_2018** (+ `prad_tcga` for Gleason) | [TCGA data dictionary](docs/tcga_prad_data_dictionary.md); see [TCGA-PRAD](#tcga-prad) |
| Cambridge / CamCap | `cambridge-gse70768` | GEO **GSE70768 / GPL10558**, Europe PMC **PMC4588396** | [Cambridge guide](docs/cambridge.md) |
| Stockholm | `stockholm-gse70769` | GEO **GSE70769 / GPL10558**, Europe PMC **PMC4588396** | [Stockholm guide](docs/stockholm.md) |
| Mortensen / Aarhus | `mortensen-gse46602` | GEO **GSE46602 / GPL570** | [Mortensen guide](docs/mortensen.md) |
| Long FFPE RNA-seq | `long-gse54460` | GEO **GSE54460**, supplementary FPKM matrix | [Long guide](docs/long.md) |
| CIT / Kamoun | `cit-e-mtab-6128` | BioStudies **E-MTAB-6128**; raw CEL + SDRF | [CIT guide and R setup](docs/cit.md) |

Supported cohorts are listed above. See the [dataset compatibility assessment](docs/dataset_compatibility.md)
for remaining cohorts and source limitations. The sections below describe Taylor unless they say otherwise.

### Sboner / GSE16560 (GPL5474)

```bash
uv run --locked sboner-gse16560                      # snapshot sboner-v2
uv run --locked sboner-gse16560 --offline            # rebuild from cached, checksummed inputs
```

Same output layout and reproducibility rules as Taylor (`data/raw|interim|processed|anndata/<snapshot>/`),
reusing the shared source cache, CSV writer, AnnData export and output protection.
Differences: GPL5474 supplies Entrez ID/symbol/RefSeq itself, so annotation is validated
rather than looked up in MyGene (6,144 probes, all retained as separate features, 6,100
unique genes); clinical data come only from GEO sample characteristics. All 281 samples
(one FFPE primary tumor per patient) are kept. Endpoints follow the shared contract (`<e>_event`, `_months`, `_eligible`, `_eligibility_reason`,
`_time_unit`): primary `pcsm` (prostate-cancer-specific mortality, encoded from the deposited
extreme-groups label: 165 lethal events, 116 indolent) and `os` (all-cause death, 206 events), both
with `followup_months`. Missing outcomes stay missing and ineligible; age, Gleason, diagnosis year,
percent cancer and fusion are parsed strictly (bad values stop the build; Gleason primary + secondary
must equal the deposited score). One sample, GSM416492, has minor Gleason pattern 6 (3+6=9); it is
retained as deposited and listed in `validation.json`. **Caveats:** lethal/indolent
proportions are set by the sampling design, not cohort incidence; the time origin is not
stated in GEO; the cohort overlaps the Swedish WW population in Setlur (GSE8402); the DASL
panel is targeted (~6k genes). Code: `src/ds_bu_bioinformatics/sboner_pipeline.py`; tests:
`tests/test_sboner_pipeline.py`.

## Quick start

Requirements: Git, Python 3.12+, and [uv](https://docs.astral.sh/uv/).
The lockfile fixes dependency versions; `.python-version` selects Python 3.12.

```bash
git clone --branch add-datasets https://github.com/BU-Spark/ds-bu-bioinformatics.git
cd ds-bu-bioinformatics
uv sync --locked
uv run --locked taylor-gse21034
```

The command downloads GEO expression/platform records, cBioPortal clinical data
and MyGene mappings, validates them, and writes the outputs below. First-time
acquisition needs network access and depends on upstream service availability.
Failed or invalid downloads retry up to three times and are not cached as valid
inputs. No manual source download is required when the services are available.

For Cambridge, after installing dependencies, run:

```bash
uv run --locked cambridge-gse70768
```

This creates the `cambridge-v2` snapshot in separate folders under the same data
roots. It keeps all 112 primary tumors, with BCR outcomes reconciled against the
original publication clinical table (19 events, 93 censored). See the
[Cambridge guide](docs/cambridge.md) before analysis.

For Stockholm, run:

```bash
uv run --locked stockholm-gse70769
```

This creates the separate `stockholm-v1` snapshot: 94 primary tumors, with 92
eligible for BCR analysis (45 events, 47 censored). Two samples with missing
clinical outcomes remain in expression and H5AD; see the [Stockholm guide](docs/stockholm.md).

For Mortensen and Long:

```bash
uv run --locked mortensen-gse46602
uv run --locked long-gse54460
```

Mortensen retains 36 tumors with 22 BCR events. Long retains 100 patients after
outcome-independent replicate selection: 49 source-coded recurrences, of which
40 are survival-eligible; nine zero-month events remain flagged for review.
Both export separate features and CSV/H5AD with full selection/error audits.
See the [Mortensen](docs/mortensen.md) and [Long](docs/long.md) guides, especially
Long's unresolved clinical-time discrepancies and zero-month events.

The commands and policies below describe Taylor unless stated otherwise.

## Output layout

```text
data/
  raw/<snapshot>/          Original source files/responses and manifest.json
  interim/<snapshot>/      feature_expression.csv.gz before mapping selection
  processed/<snapshot>/
    expression.csv.gz      Original feature IDs × selected sample IDs
    genes.csv              Feature ID, gene ID/symbol and RefSeq annotations
    samples.csv            Sample, GEO and cBioPortal patient identifiers
    clinical.csv           Clinical covariates and nullable DFS outcomes
    feature_mapping.csv.gz Every feature/candidate and its mapping decision
    sample_audit.csv       Every GEO profile and its selection reason
    metadata.json          Source/code hashes, schemas and processing policies
    validation.json        Counts, missingness, alignment and round-trip checks
    checksums.sha256       SHA-256 hashes of the final CSV/metadata files
  anndata/<snapshot>/
    expression.h5ad        Samples × features, clinical data and gene annotations
    metadata.json          Same source/processing metadata as the CSV export
    checksums.sha256       SHA-256 hashes of H5AD and metadata
```

The default snapshot is `taylor-v2`. Data, caches and generated outputs are
ignored by Git. The synthetic files in `tests/fixtures/` are committed.

In H5AD, `X` is the transposed expression matrix, `obs` contains clinical/sample
fields indexed by sample ID, and `var` contains gene/RefSeq annotations indexed
by original feature ID. Gene IDs and symbols may repeat; feature IDs stay unique.
`uns['metadata_json']` stores source and processing provenance as JSON. CSV and
H5AD reloads must preserve expression values, labels and annotation alignment
before outputs are published. See the [data dictionary](docs/data_dictionary.md)
for field meanings and clinical limitations.

## Reproduce a snapshot

Every cached input has its request, retrieval timestamp, byte size and SHA-256
recorded in `data/raw/<snapshot>/manifest.json`. Repeated runs verify and reuse
these inputs. Existing differing outputs are protected rather than overwritten.
After acquisition, rebuild without network access:

```bash
uv run --offline --locked taylor-gse21034 --offline --snapshot taylor-v2
```

`uv --offline` prevents dependency downloads; the pipeline's `--offline` prevents
HTTP acquisition and fails on missing or corrupt inputs. For an independent
rebuild, use **all three** alternate output roots:

```bash
uv run --offline --locked taylor-gse21034 --offline --snapshot taylor-v2 \
  --output-root /tmp/taylor-check/processed \
  --interim-root /tmp/taylor-check/interim \
  --anndata-root /tmp/taylor-check/anndata
```

To reproduce the same biological inputs on another machine, transfer the complete
`data/raw/taylor-v2/` directory, including its manifest, then run offline with the
locked environment. Also install dependencies before disconnecting. The
[reference manifest](docs/provenance/taylor-v2-manifest.json) records the 58 inputs
used for the verified cohort; it contains requests and checksums, **not the input
payloads**. A fresh online download can differ if upstream annotations or clinical
records change, even under the same snapshot name in a new checkout. See
[provenance notes](docs/provenance/README.md).

To refresh sources intentionally, use a new input snapshot:

```bash
uv run --locked taylor-gse21034 --snapshot taylor-v3
```

If pipeline code changes and existing output metadata differs, keep the raw
snapshot and choose new output roots. CSV/gzip outputs are deterministic;
checksum manifests allow file-level comparisons. H5AD also has semantic reload
checks. Code hashes in metadata intentionally change when implementation changes.
Taylor records the current package code hashes on a fresh export. When rerunning
with both published CSV and AnnData exports present, Taylor uses their original
code provenance for verification and still requires every rebuilt file to match.

## Interpretation and policies

- **Expression:** retain finite deposited values. GEO reports RMA background
  adjustment and quantile normalization, but the exact deposited log base remains
  unverified. No additional normalization, log transform, imputation, expression
  filtering or batch correction is applied. No raw CEL quality control is done.
- **Mapping:** GPL10264 maps features to RefSeq. Preserve original accessions;
  strip only a terminal numeric version for lookup. Cache complete MyGene human
  responses. Require exact RefSeq RNA agreement and a consistent Entrez/symbol
  identity. Exclude ambiguous, unmapped, invalid, truncated or conflicting
  mappings with reasons. Never select by highest score. Symbols are MyGene
  annotations, not claimed HGNC-approved labels.
- **Separate features:** retain each accepted feature and its deposited values,
  including features sharing gene IDs or symbols. No averaging is applied.
- **Samples:** require GEO `Primary tumor` and cBioPortal `Primary` when supplied.
  Match exact sample IDs and verify patient associations. Missing joins,
  conflicting records and multiple selected primary profiles per patient fail.
  Identical duplicate clinical records are collapsed and counted.
- **Clinical:** retain original covariates with `patient__`/`sample__` prefixes.
  Missing values remain missing. Source definitions and provenance are recorded.
- **DFS:** map source `0:DiseaseFree` to event 0 and `1:Recurred` or
  `1:Recurred/Progressed` to event 1; unknown codes fail. Times remain in months.
  Missing outcomes remain nullable. Negative/non-finite times fail; zero times
  are retained but flagged ineligible. The exact recurrence/death criteria and
  surgery-specific time origin remain unverified. DFS is not relabeled BCR.

The verified frozen `taylor-v2` build contains **131 primary profiles, 25,900
separate mapped features, 19,306 unique genes, 27 DFS events and 104 disease-free
observations**. Counts are calculated, not imposed. `gene_rows` in the validation
report counts unique Entrez genes; `retained_feature_rows` counts expression rows.
Legacy counts are comparisons, not selection targets. Structural validity does
not resolve the remaining expression-scale or DFS-definition limitations.

## TCGA-PRAD

```bash
uv run --locked tcga-prad                                   # snapshot tcga-prad-v2
uv run --locked tcga-prad --module-file modules.csv         # also writes module_coverage.csv
uv run --offline --locked tcga-prad --offline --snapshot tcga-prad-v2   # rebuild from cache
```

Inputs are the cBioPortal datahub RSEM matrix (`data_mrna_seq_v2_rsem.txt`, ~77 MB, pinned to
datahub commit `86690e1`),
the PATIENT/SAMPLE clinical records and sample inventory from the cBioPortal API, and
Gleason fields from the Firehose Legacy study. All of them are cached and checksummed
the same way as Taylor's. Outputs go to `data/processed/<snapshot>/` and
`data/anndata/<snapshot>/`, using the Taylor file names. Two names differ:

- `feature_audit.csv` replaces `feature_mapping.csv.gz`, because there is no RefSeq mapping step.
- `qc_summary.md` is added as a readable summary.

TCGA-specific policies (CampLab client meetings and the 11/06 team meeting):

- **Expression:** values are **log2(RSEM + 1)**. PanCanAtlas batch correction leaves
  small negative values, and these are clipped to 0 first. The unclipped RSEM values are
  in the h5ad layer `rsem`. Rows with a missing value in any selected sample
  (batch-correction gaps) are removed, not imputed. No rows are collapsed: `feature_id` is the Entrez ID, or
  `<entrez>_row<N>` when an ID repeats in the source.
- **Samples:** barcode characters 1–12 are the patient and 14–15 the sample type. The
  pipeline:
  1. keeps type `01` (primary tumor);
  2. keeps one sample per patient (the lowest barcode, which is outcome-independent);
  3. excludes expression QC outliers (|z| > 3 on median expression or genes detected).

  The barcode patient must match the cBioPortal inventory, and cBioPortal `SAMPLE_TYPE`
  must be `Primary`.
- **Endpoints:** **PFI** (`pfi_*`, from `PFS_STATUS`/`PFS_MONTHS`) is primary. OS, DSS
  and DFI follow the same nullable event/time/eligibility pattern as Taylor `dfs_*`, and
  unknown status codes fail. PFI is not relabeled BCR or DFS. Event counts in
  `validation.json` include eligible records only.
- **Covariates:** `age` and `gleason_score` fail on non-numeric or out-of-range values,
  the same as endpoints; only missing-value tokens become blank.

The November 11, 2026 `tcga-prad-v2` build contains **484 primary tumors (one per
patient; 9 QC outliers excluded), 19,088 feature rows, 19,064 unique genes and 91
PFI events** (median observed time 25.8 months). Offline rebuilds are byte-identical,
including the h5ad. The [reference manifest](docs/provenance/tcga-prad-v2-manifest.json)
records the inputs. `tcga-prad-v2` replaces `tcga-prad-v1`: the inputs and data files are
identical, but v1's cached RSEM request used the unpinned datahub `master` URL, so the
current code cannot rebuild it.

## Checking all datasets

Each dataset has its own pipeline, because source formats, gene identifiers, sample
rules and endpoints differ (see [dataset compatibility](docs/dataset_compatibility.md)).
What they share is the output format. After running the pipelines, check every
snapshot in `data/processed/` and `data/anndata/` against that format:

```bash
uv run --locked check-outputs                          # every snapshot
uv run --locked check-outputs --snapshot taylor-v2     # one snapshot (repeatable)
uv run --locked check-outputs --json report.json       # also save full results
```

It does not rerun any pipeline. For each snapshot it checks:

- required files, plus a feature audit (`feature_mapping.csv.gz` or `feature_audit.csv`);
- `checksums.sha256` in both folders;
- `expression.csv.gz` is feature × sample with unique labels and finite values;
- sample order is the same in `samples.csv`, `clinical.csv` and the selected rows of
  `sample_audit.csv`, with one sample per patient;
- feature order is the same in `genes.csv`;
- every `<endpoint>_event` has `_months`, `_eligible`, `_eligibility_reason` and
  `_time_unit` columns, with valid values; all time units must be `months`, and
  non-missing follow-up times must be finite;
- `validation.json`, `metadata.json` and `expression.h5ad` agree with the CSVs.

It prints one row per snapshot, with counts computed from the files rather than
from `validation.json`, whose count names differ by dataset. It exits with status 1
if any snapshot fails. Warnings flag differences that don't break the format but
matter for combined analysis.

Result on November 13, 2026 (Cambridge and Stockholm built from the `add-datasets` branch):

| Dataset | Status | Samples | Feature rows | Unique genes | Expression min / median / max | Endpoints (eligible, events) |
|---|---|---:|---:|---:|---|---|
| taylor-v2 | PASS | 131 | 25,900 | 19,306 | 11.72 / 224.07 / 30762.05 | DFS 131, 27 |
| tcga-prad-v2 | PASS | 484 | 19,088 | 19,064 | 0.00 / 7.95 / 20.62 | PFI 484, 91; OS 484, 10; DSS 482, 5; DFI 330, 30 |
| cambridge-v2 | PASS | 112 | 41,610 | 29,604 | 4.94 / 6.63 / 15.38 | BCR 112, 19 |
| stockholm-v1 | PASS | 94 | 41,615 | 29,608 | 5.30 / 6.12 / 14.50 | BCR 92, 45 |

Warning: Taylor's deposited values are not on a log scale, while the other three are
log2. Decide how to transform them before PCA module scoring. Endpoints also differ
(DFS, PFI, BCR) and are not relabeled, so compare module scores within each
dataset's own endpoint.

## Tests and report

```bash
uv run --locked python -m unittest discover -s tests -v
uv run --locked jupyter notebook notebooks/taylor_gse21034_report.ipynb
uv run --locked jupyter notebook notebooks/tcga_prad_report.ipynb   # TCGA-PRAD
```

The tests use synthetic inputs and mocked network responses. They cover sample
identity/order, mapping ambiguity, separate features with repeated gene labels,
nullable clinical outcomes, malformed inputs, cache corruption, pagination,
interrupted acquisition, output protection and CSV/H5AD round trips. GitHub
Actions runs the tests with the lockfile and Python 3.12; no live cohort download
is required for tests.

The notebook reads local outputs, verifies checksums/alignment, and demonstrates
separate features sharing a gene. Run the pipeline first. To execute the report
without opening Jupyter:

```bash
uv run --locked jupyter nbconvert --to notebook --execute \
  notebooks/taylor_gse21034_report.ipynb --output /tmp/taylor-report.ipynb
```

## Repository layout

- `src/ds_bu_bioinformatics/`: CLI, acquisition, parsing, annotation, clinical
  matching, validation and AnnData export. `tcga_pipeline.py` holds the TCGA-PRAD
  pipeline and reuses the shared cache, clinical paging, CSV and AnnData writers.
- `tests/`: synthetic fixtures and automated tests.
- `notebooks/`: local inspection report; generated cell outputs are not committed.
- `docs/`: data dictionary, cohort compatibility assessment and input provenance.
- `pyproject.toml`, `uv.lock`: package/CLI definition and locked dependencies.

Superseded project templates were replaced by this workflow and documentation.
Previous averaged outputs remain local under `data/processed-averaged/` and are
not inputs to the feature-preserving workflow. License: [GPL-3.0](LICENSE).
