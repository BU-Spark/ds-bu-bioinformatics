# Mortensen / Aarhus GSE46602

## Run and outputs

```bash
uv sync --locked
uv run --locked mortensen-gse46602
```

Default snapshot: `mortensen-v1`. The pipeline automatically downloads the processed
GEO Series Matrix and GPL570 annotation. It writes the standard expression,
genes, samples, clinical, sample audit, feature mapping, metadata, validation and
checksum files under `data/processed/mortensen-v1/`; selected source expression
under `data/interim/mortensen-v1/`; and `expression.h5ad`, metadata and checksums
under the separate `data/anndata/mortensen-v1/` folder. Raw inputs and their request
manifest remain in `data/raw/mortensen-v1/`. Generated data are ignored by Git.

H5AD `X` is samples × probes, `obs` contains clinical/sample fields, `var` contains
probe annotations, and `uns['metadata_json']` records provenance and policies.

## Policies

- Select `tissue: prostate tumor` and `Prostate cancer cells, PC*` titles.
  Exclude both ordinary benign glands and the explicitly labeled glands adjacent
  to cancer. Unknown tissue/title combinations fail for review.
- Sample IDs are GEO GSM accessions. Patient IDs are the `PC*` labels in the GEO
  titles, not independently linked registry IDs. Duplicate selected subjects fail.
- Preserve deposited finite RMA-normalized values. No additional transformation,
  normalization, filtering by expression magnitude, imputation or probe averaging.
- Use deposited GPL570 human Entrez IDs, symbols and RefSeq annotations. Require
  one gene ID/symbol identity; exclude multi-gene probes, invalid annotations and
  inconsistent gene/symbol identities. Keep distinct probes mapping to the same
  accepted gene separate. This is frozen platform annotation, not a fresh MyGene
  or HGNC revalidation.
- If expression is missing in any selected primary, exclude that whole probe and
  record the decision. Missing values only in excluded benign samples do not
  affect primary-tumor exports. Invalid numeric text/infinities fail.
- Map source `bcr`: `YES → 1`, `NO → 0`; use deposited `bcr_free_time` in months
  for event or censoring duration. The primary paper describes postoperative
  follow-up in months and biochemical failure as PSA >0.2 ng/ml. No time conversion
  or independent clinical adjudication is performed.
- Retain missing outcomes and zero times with explicit eligibility reasons;
  never assume missing means censored. Negative/nonfinite clinical times fail.

The frozen November 13, 2026 build contains **36 primary tumors, 22 eligible BCR
events, 14 censored observations, 41,692 separate probes and 20,462 unique Entrez
IDs**. All 36 primaries are BCR-eligible. The 14 benign profiles remain in the
sample audit. The full source matrix has 54,675 probes; annotation policy retains
41,692, with no missing-expression exclusions in the selected tumors. Counts are
computed, never enforced as selection targets.

## Reproduce and verify

```bash
uv run --offline --locked mortensen-gse46602 --offline \
  --output-root /tmp/mortensen-check/processed \
  --interim-root /tmp/mortensen-check/interim \
  --anndata-root /tmp/mortensen-check/anndata
uv run --locked python -m unittest discover -s tests -v
```

For another machine, transfer the complete raw snapshot and install locked
dependencies before disconnecting. The reference manifest in
`docs/provenance/mortensen-v1-manifest.json` contains request details and hashes,
not payloads. Existing differing outputs are protected across all formats. Use
new output roots after code changes or a new snapshot for intentional input refreshes.

The CSV reload checks expression values and all exported tables. The H5AD reload
checks numeric values, row/column order, annotations and metadata. Tests include
ambiguous annotations, repeated gene IDs, missing outcomes, benign exclusion,
missing probes, malformed values and offline reproduction.

Sources: [GSE46602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46602),
[GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570),
[primary publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC4629186/).
