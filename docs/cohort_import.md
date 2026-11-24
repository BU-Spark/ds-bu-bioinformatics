# Optional curated local-file import

All supported cohorts also have dedicated commands; see
[pipeline coverage](pipeline_coverage.md). This utility is an alternative for
explicitly curated inputs and is not counted as a dataset-specific pipeline.

`import-cohort` provides the common CSV/AnnData pipeline for research files that have
already been acquired and **curated into an explicit schema**. It is not an automated
parser for each original study, a source-data discovery tool, or proof of survival
readiness. Missing clinical follow-up and controlled-access authorization remain real
prerequisites. The implementation has synthetic integration coverage, not full-source
validation for these cohorts.

Supported identities:

| `--cohort` | Required accession |
|---|---|
| `nakagawa` | GSE10645 |
| `cpc-gene` | GSE107299 |
| `cancermap` | GSE94767 |
| `erho` | GSE46691 |
| `karnes` | GSE62116 |
| `den` | GSE72291 |
| `klein` | GSE62667 |
| `moscow` | GSE229904 |
| `dkfz-eopc` | EGAD00001004791 |
| `cpgea` | PRJCA001124 |

Do not substitute another cBioPortal/DKFZ cohort. Do not assume literature cohort
sizes equal expression sample counts. For controlled files, acquire them through the
appropriate research access process before preparing a bundle; this importer only
reads explicitly supplied local files.

## Bundle files

Put these files next to a JSON manifest, using relative paths within its directory:

| Logical input | Required schema |
|---|---|
| `expression` | Gzipped CSV: first column `feature_id`, then unique source sample IDs. Preserve original source feature rows and numeric values. Include every source assay listed in the sample audit. |
| `genes` | CSV: unique `feature_id`, `gene_id`, `gene_symbol`. Exact coverage of expression rows. Canonical IDs can be blank when unresolved; do not invent Entrez IDs. Ambiguous/conflicting gene annotations are excluded and audited. |
| `samples` | CSV: unique nonempty `sample_id`, `patient_id`, `platform_id`, one selected assay per patient. Preserve source identity; describe repeated-region/panel selection explicitly. |
| `clinical` | CSV: exact selected `sample_id` coverage and `<endpoint>_event`, `<endpoint>_time`, `<endpoint>_unit` for each declared endpoint. Events are `0`/`1` or missing; times are individual durations or blank; units are `months`/`years` even if time is missing. Extra source covariates can be retained. |
| `sample_audit` | CSV: unique source `sample_id`, `selected` true/false, and nonempty `selection_reason`. Must cover every expression column; selected rows must follow `samples.csv` order. |

Clinical rows may arrive in a different order and are aligned by an exact key.
Gene rows are likewise aligned by exact feature ID. Patient duplication, missing or
extra clinical IDs, a partial sample audit and mismatched annotations fail.

These inputs must be derived from original research files using a reviewed, documented
curation procedure. An audit table records supplied selection; software cannot establish
whether a manually supplied patient crosswalk or clinical definition is true. For
Nakagawa, prefer its dedicated paired-panel adapter; a generic bundle would require
an independently curated combined matrix and assay crosswalk.

## Manifest

Start with [the manifest example](templates/cohort_bundle.example.json). Replace every
placeholder and supply each actual file's SHA-256 (`shasum -a 256 path/to/file`). The
manifest requires:

- `schema_version: 1`, exact `study_accession`, `organism: "Homo sapiens"`, and `platforms`;
- `expression_kind`, plus a precise human-readable `expression_scale`;
- `source_description`, `annotation_source`, and `selection_policy` describing original
  sources, processing, patient mapping and selection independently of outcomes;
- `primary_endpoint`, and an `endpoints` object documenting each endpoint's definition,
  time origin, censoring and sample linkage;
- five `files` entries with relative `path`, actual `sha256` and original `source` citation.

Accepted expression kinds:

| Kind | Export | Module-scoring transform |
|---|---|---|
| `log_microarray` | Original supplied values | None |
| `linear_microarray` | Original supplied values | log2(x + 1) |
| `tpm` / `fpkm` | Original supplied values | log2(x + 1) |
| `raw_counts` | Original supplied integer counts | log2(CPM + 1), using full selected-source library totals frozen before feature filtering |

Raw counts must be complete, finite, nonnegative integers with positive library totals.
This basic log-CPM option is explicit and reproducible; it is not a claim that library
size scaling resolves every compositional/batch effect. Counts from miRNA assays must
not be mixed with mRNA. Other count normalizations require separately documented work.

Missing times remain missing and ineligible. A structural validation PASS means the
export is internally consistent, not that it is suitable for Cox modeling. Modules
can still be scored without fitting a survival model. Case-control sampling and
cohort overlap require analysis-specific treatment beyond this importer.

## First import and offline rebuild

```bash
uv run --locked import-cohort --cohort erho --bundle /path/to/bundle/manifest.json
uv run --locked check-outputs --snapshot erho-import-v1
uv run --locked module-scores --snapshot erho-import-v1
```

Every input, including the manifest, is copied into the checksummed raw snapshot.
Changed files, mismatched hashes and path traversal fail. A cached replay does not
need the original bundle or a network connection:

```bash
uv run --locked import-cohort --cohort erho --offline \
  --output-root data/rebuilt/processed --interim-root data/rebuilt/interim \
  --anndata-root data/rebuilt/anndata
```

Use a new snapshot when changing inputs. The optional `--clinical-file` and
`--clinical-provenance` flags support [additional verified clinical linkage](clinical_linkage.md)
without overwriting original source fields.
