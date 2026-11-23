# Additional public GEO pipelines

Each cohort below now has its own native acquisition command. These commands fetch
original GEO metadata, processed expression and annotation inputs, perform the cohort's
patient/specimen selection, and generate aligned CSV/AnnData exports. They do not
require users to prepare the generic five-file import bundle.

**Validation status:** synthetic source-format integration and offline-rebuild tests
cover all seven workflows. Complete real-source acquisition, platform/source-schema
coverage and cohort counts have not yet been validated in this execution environment.
A successful fixture test is not a claim that a real-source run has completed.

| Cohort | Command | Snapshot | Native inputs and selection |
|---|---|---|---|
| CancerMap | `cancermap-gse94767` | `cancermap-v1` | GSE94767 Series Matrix, GPL5175 annotation; human malignant tissue, explicit source donor labels, one region per donor |
| CPC-GENE | `cpc-gene-gse107299` | `cpc-gene-v1` | Both platform-specific GSE107299 metadata matrices plus deposited processed TSV; exact CPCG/GSM patient mapping |
| Erho / Mayo | `erho-gse46691` | `erho-v1` | GSE46691 metadata/matrix, deposited quantile-normalized supplementary PSR matrix when the Series Matrix is empty, GPL5188 |
| Karnes / Mayo | `karnes-gse62116` | `karnes-v1` | GSE62116 processed Series Matrix and GPL5188 |
| Den / Thomas Jefferson | `den-gse72291` | `den-v1` | GSE72291 processed Series Matrix and GPL5188; retain postoperative-RT context |
| Klein / Cleveland | `klein-gse62667` | `klein-v1` | GSE62667 processed Series Matrix and GPL5188; retain source SCAN processing |
| Moscow / EMB-RAS | `moscow-gse229904` | `moscow-v1` | GSE229904 metadata plus deposited mRNA count TSV; tumor mRNA assays only, exclude normals/miRNA |

## Running and rebuilding

```bash
uv sync --locked
uv run --locked cancermap-gse94767
uv run --locked cpc-gene-gse107299
uv run --locked erho-gse46691
uv run --locked karnes-gse62116
uv run --locked den-gse72291
uv run --locked klein-gse62667
uv run --locked moscow-gse229904
```

Run each command once with network access to freeze every input. For example:

```bash
uv run --locked check-outputs --snapshot moscow-v1
uv run --locked moscow-gse229904 --offline \
  --output-root data/rebuilt/processed --interim-root data/rebuilt/interim \
  --anndata-root data/rebuilt/anndata
```

Compare original/rebuilt file hashes before recording a validated real cohort.
Use a new snapshot when input requests or source content change. Existing outputs are
protected from replacement by incompatible results. New code files change provenance
hashes: rebuilding an old export after a code change may require separate output roots.

The Series Matrix/TSV readers convert numeric data in batches rather than constructing
one Python float object per matrix element. Exon-array matrices can nevertheless be
very large. Budget substantial memory, disk space and runtime for full feature-level
CSV/H5AD exports, particularly Erho. No automatic gene averaging is used to reduce size.

## Shared annotation and missing-data rules

Affymetrix workflows retain original feature IDs and full source gene-assignment text.
A source assignment must name one gene symbol and provide a valid RefSeq anchor;
cached human mapping responses must return that exact RefSeq. Ambiguous, unmapped
and conflicting identities are quarantined, with mapping candidates exported.
GPL5188 PSRs are never collapsed to transcript clusters. GPL5175 source gene-target
summaries remain separate; raw probes already summarized by the source cannot be
recovered from the processed matrix.

CPC-GENE uses its deposited Entrez/Brainarray identifiers; Moscow detects a consistent
Entrez, Ensembl or symbol namespace. Exact human mappings are cached and replayed
offline. Source identifiers and distinct rows remain in the audit even when their
canonical annotations repeat. A row with missing expression in any selected subject
is excluded and audited; there is no zero fill or expression-dependent patient selection.

## Cohort-specific issues

### CancerMap

Select Homo sapiens malignant tissue using explicit tissue/donor characteristics.
Supported source-field aliases are listed in `cancermap.py`; unknown source fields
or tissue values fail and print the fields that need review. Never assume one patient
per region/array. Among multiple eligible regions for a donor, the lowest GSM identifier
is selected independently of outcomes and the other regions remain in the audit.
This is a reproducibility rule, not evidence that the chosen region is clinically
dominant. Region sensitivity analysis remains a downstream research decision.

### CPC-GENE

Original batch/platform fields remain in clinical data. Processed columns must map
exactly to source CPCG/GSM labels. The GEO deposit documents an upstream average of
CPCG0398 and CPCG0398_rep; the pipeline does not average again or count the technical
replicate as another patient. Deposited probe summarization and batch correction are
not reversed. Do not equate the proteomic-paper subset size with the complete expression
cohort. Clinical tables must be joined by verified source CPCG/GSM crosswalks.

### Erho, Karnes, Den and Klein

Require primary-tumor source titles, unique study-qualified subjects and GPL5188.
Erho's supplementary matrix aliases use exact source descriptions or GSMs. A processed
column that cannot be mapped uniquely fails rather than being positionally matched.
The deposited metastatic-event label is retained when present. Individual durations
remain missing until verified. Den also declares BCR as a secondary endpoint and
preserves the postoperative-radiotherapy context. Outcomes and times are not fabricated
from Gleason score, PSA values, treatment labels or published event totals.

These Mayo/Decipher studies can overlap and include selected sampling designs.
Study-qualified IDs keep separate exports distinct; they do not establish independence
or justify pooling cohorts without a patient-overlap/design review.

### Moscow

Only `*T RNA` primary tumor samples enter expression/PCA. Normal `*N RNA` and miRNA
assays are audited and excluded; exact GSM, full title or specimen-core IDs link count
columns to RNA metadata. Primary RNA profiles absent from the submitted count matrix
are explicitly counted in the audit, not silently assumed to have expression.

Source `bcr`/`brf` codes map to known recurrence classes. Individual follow-up remains
missing. Raw counts must be complete nonnegative integers with a positive library total.
Exports retain counts unchanged and freeze library totals from the full selected-source
matrix **before** annotation filtering. Downstream normalization and analysis are outside
this dataset-pipeline change; raw counts should not be treated as normalized expression.

## Clinical supplementation and analysis

Every new public command accepts `--clinical-file` and `--clinical-provenance` through
[the exact-linkage interface](clinical_linkage.md). It fills verified missing endpoints
without overwriting contradictory source values. Missing follow-up does not block the
expression export. Consumers must separately verify any subsequent analysis eligibility.
A structural `check-outputs` PASS is not a declaration of survival readiness.

Source studies:
[GSE94767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94767),
[GSE107299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107299),
[GSE46691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46691),
[GSE62116](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62116),
[GSE72291](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72291),
[GSE62667](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62667), and
[GSE229904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229904).
