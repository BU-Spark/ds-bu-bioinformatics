# Long FFPE RNA-seq GSE54460

## Run and output layout

```bash
uv sync --locked
uv run --locked long-gse54460
```

Default snapshot: `long-v1`. This downloads two processed source files from GEO:

1. `GSE54460_series_matrix.txt.gz`: clinical metadata for 106 assays; no expression
   rows. Repeated processing headers are preserved by the shared GEO reader.
2. `GSE54460_FPKM.txt.gz`: original supplementary FPKM matrix, including its
   embedded sample/clinical metadata before the Entrez/SYMBOL table.

The standard expression/genes/samples/clinical/feature-mapping/sample-audit,
metadata, validation and checksum files are under `data/processed/long-v1/`.
`data/interim/long-v1/feature_expression.csv.gz` retains all source feature rows
for selected assays. H5AD and its metadata/checksums are separate under
`data/anndata/long-v1/`. Raw source payloads and the download manifest remain
under `data/raw/long-v1/`. Generated files are ignored by Git.

Two additional audit files are included:

| File | Contents |
|---|---|
| `clinical_source_comparison.csv` | Every assay's GEO and embedded-FPKM status/time, linkage method, selected flag and numerical agreement flags |
| `expression_value_audit.csv` | Every recognized missing/error cell: original feature-row ID, expression assay ID, source token, reason, GEO accession and whether the assay was selected |

H5AD uses samples × source-feature rows in `X`, clinical fields in `obs`, original
source annotations in `var`, and JSON provenance in `uns`. No assay or feature is
averaged.

## Sample identities and replicates

The source study consists of FFPE radical-prostatectomy tumor material. Require
GEO tissue `prostate` and the expected `Prostate_<external ID>` title (including
explicit `-2` replicate suffix). The adapter does not claim GEO supplies a separate
Primary sample-type label; the study design provides that cohort interpretation.

The 106 expression assay IDs differ from GEO `cgsr id` in 12 cases:

- Match 94 by exact cgsr ID and verify the external patient ID independently.
- Match the 12 aliases to the unique **remaining** GEO assay with the same external
  ID. Require a complete one-to-one assay bijection. Ambiguity, duplicate IDs or
  conflicting external IDs fail; no matches are inferred from column position.

Keep one assay per external patient ID: choose the greatest GEO uniquely mapped
read count, with the lowest GSM accession as tie-break. This is independent of
outcomes and expression values. Invalid/missing quality counts fail rather than
silently selecting another assay. Six patients have two assays; six lower-priority
assays are audited and excluded. Original selected sample order follows GEO.
Conflicting known BCR statuses within a patient fail.

## Expression rows and annotations

Retain deposited **FPKM**, with no extra normalization or log transform. FPKM is
not a read-count matrix or directly interchangeable with log-scale microarray
signals. Negative values, infinities and unknown nonnumeric tokens fail.

The source contains 23,281 rows but only 23,207 distinct Entrez IDs. A total of 139
rows share an Entrez ID with another row; some have different deposited labels.
Preserve them all separately. Feature keys are `GSE54460:line:<N>`, where N is the
original one-based physical file line. `gene_id`, `gene_symbol`, `source_row` and
`gene_id_repeated_in_source` remain annotations. These legacy gene assignments
are not claimed current, internally reconciled or HGNC-approved; no fresh
mapping or gene-symbol substitution is performed.

One source cell is `#NUM!`: MIR4426 / Entrez 100616345 in assay `CM.4-0021`.
It is recorded in `expression_value_audit.csv`, and its entire feature row is
excluded because that assay is selected. Recognized missing/error cells are never
filled with zero. A missing value only in a discarded replicate does not remove
a feature complete across selected assays. Row IDs are assigned before filtering
so they do not change when another feature is excluded.

The verified build retains **23,280 separate source rows, 23,206 unique deposited
Entrez IDs and 100 patients**, one selected assay per patient.

## Clinical source and remaining source limitations

Use **GEO sample characteristics consistently** for source BCR status, event time
and censoring follow-up. Keep all original `geo__` fields and embedded `fpkm__`
fields side by side. This is an explicit reproducible source choice, not proof
that GEO is the most up-to-date clinical authority. The primary paper and source
deposits require further clinical adjudication to establish an authoritative
resolution of the differences; the pipeline does not mix or silently repair them.

All 106 matched assay BCR statuses agree between GEO and the FPKM table. Their
clinical times differ in 40 assay event-time entries and 62 follow-up entries.
Among selected assays, 36 event-time and 62 follow-up entries differ. Some
variation may reflect precision or updates; its cause remains unverified.
Comparison uses absolute tolerance 1e-8 months (no relative tolerance), with
missing versus missing considered agreement. Original text is also preserved.

Map GEO `bcr`: `1 → event`, `0 → non-event`. Use `months to bcr` for events and
`months total f/u` for censoring. Do not use the FPKM embedded time as fallback.
Missing clinical values stay nullable; negative/nonfinite/malformed times fail.
An event later than recorded follow-up is flagged ineligible.

Of the 100 selected patients, **49 have source-coded recurrence and 51 are
non-events**. Nine recurrence patients have zero-month event times. The zeroes
also occur in the embedded FPKM source; no positive event time can be inferred.
Keep all nine patients and source events but flag `zero_time_requires_review`.
No pseudotime is added. This leaves **91 survival-eligible patients: 40 BCR events
and 51 censored**. `validation.json` clearly distinguishes eligible `bcr_events`
from `source_coded_bcr_events`. Counts are observations, not forced targets.

## Reproduce and test

```bash
uv run --offline --locked long-gse54460 --offline \
  --output-root /tmp/long-check/processed \
  --interim-root /tmp/long-check/interim \
  --anndata-root /tmp/long-check/anndata
uv run --locked python -m unittest discover -s tests -v
```

Transfer the entire `data/raw/long-v1/` cache and install locked dependencies for
an offline rebuild on another machine. `docs/provenance/long-v1-manifest.json`
records source requests/hashes, not the payloads. Choose new output roots after
implementation changes. Existing differing outputs are protected across formats.
CSV reload checks all tables and exact numeric expression; H5AD reload checks
values, clinical/feature annotations, sample order and provenance.

Tests cover assay aliases, ambiguous/conflicting identity, quality ties, shuffled
source order, zero-time retention, clinical-source disagreement, duplicate gene
IDs, the numeric-error cell, missing expression only in discarded assays,
malformed source values, empty error audits and byte-identical offline rebuilds.

Sources: [GSE54460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54460),
[primary publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC4058362/).
