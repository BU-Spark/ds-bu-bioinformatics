# Taylor/MSKCC prostate cancer pipeline

Reproducible acquisition, annotation and clinical matching for **GSE21034 /
GPL10264** and cBioPortal **prad_mskcc**. The pipeline selects primary tumors,
checks sample-to-patient identities, preserves separate features sharing a gene
mapping, and exports aligned CSV and AnnData (`.h5ad`) datasets.

**Supported dataset: Taylor/MSKCC only.** Other cohorts require source,
annotation and clinical adapters; see the preliminary
[dataset compatibility assessment](docs/dataset_compatibility.md).

## Quick start

Requirements: Git, Python 3.12+, and [uv](https://docs.astral.sh/uv/).
The lockfile fixes dependency versions; `.python-version` selects Python 3.12.

```bash
git clone --branch mansur_orazkhan https://github.com/BU-Spark/ds-bu-bioinformatics.git
cd ds-bu-bioinformatics
uv sync --locked
uv run --locked taylor-gse21034
```

The command downloads GEO expression/platform records, cBioPortal clinical data
and MyGene mappings, validates them, and writes the outputs below. First-time
acquisition needs network access and depends on upstream service availability.
Failed or invalid downloads retry up to three times and are not cached as valid
inputs. No manual source download is required when the services are available.

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

## Tests and report

```bash
uv run --locked python -m unittest discover -s tests -v
uv run --locked jupyter notebook notebooks/taylor_gse21034_report.ipynb
```

The 32 tests use synthetic inputs and mocked network responses. They cover sample
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
  matching, validation and AnnData export.
- `tests/`: synthetic fixtures and automated tests.
- `notebooks/`: local inspection report; generated cell outputs are not committed.
- `docs/`: data dictionary, cohort compatibility assessment and input provenance.
- `pyproject.toml`, `uv.lock`: package/CLI definition and locked dependencies.

Superseded project templates were replaced by this workflow and documentation.
Previous averaged outputs remain local under `data/processed-averaged/` and are
not inputs to the feature-preserving workflow. License: [GPL-3.0](LICENSE).
