# CIT / Kamoun (E-MTAB-6128)

CIT is an Affymetrix Human Gene 2.0 ST cohort. BioStudies supplies raw CEL
files and the sample/clinical SDRF; it does not supply a processed matrix.
The CIT adapter runs RMA locally and then uses the same validated CSV and H5AD
export layout as the other cohorts.

## Run

The tested setup on Apple Silicon uses the pinned Intel macOS environment
under Rosetta. It does not require replacing the machine's Python or R:

```bash
micromamba create --platform osx-64 -p "$PWD/.codex-build/cit-r" \
  -f env/cit-osx-64.lock
uv sync --locked
CIT_RSCRIPT="$PWD/.codex-build/cit-r/bin/Rscript" uv run --locked cit-e-mtab-6128
```

The R script checks R 4.5.3 and the versions of oligo, platform design,
transcript-cluster annotation and relevant dependencies. Set `CIT_RSCRIPT`
to the matching executable, or expose it on PATH. The lockfile is for Intel
macOS; other platforms need a matching environment and separate verification.

The first run acquires approximately 3.7 GB of raw CELs. BioStudies inventory
and SDRF identities must agree; sizes and hashes are recorded in the raw
manifest. All 141 deposited arrays are included in the same normalization
batch. Tumor selection occurs after normalization.

## Cohort and expression policies

- The SDRF contains 141 arrays from 117 distinct patient IDs: 101 tumor arrays
  from 101 distinct patients and 40 non-tumor arrays. The paper's larger
  multi-omics cohort is not the expression deposit.
- Use the explicit `Tumoral tissue` description and keep the SDRF source,
  assay, CEL filename and patient identifiers. Normal arrays are audited but
  excluded from the final expression and clinical outputs.
- Cross-check BCR status and follow-up between the SDRF characteristics and
  factor values. Use deposited months, without conversion or imputation.
  Of the 101 tumors, 77 are survival-eligible (13 events, 64 censored).
  Twenty have unknown outcomes, and four non-events have zero-month follow-up;
  all 24 remain in expression with ineligibility reasons.
- Use `oligo::rma(target="core")` to summarize raw CEL probes into log2
  transcript-cluster intensities. This is a local RMA reprocessing, not an
  assertion of equivalence to the original paper's processed matrix.
- Keep transcript-cluster IDs as feature keys. Annotation comes from the
  pinned `hugene20sttranscriptcluster.db` package. Unmapped or ambiguous
  annotations are audited and excluded; clusters sharing a gene remain
  separate. No cross-cluster averaging, extra transformation or imputation.

The verified full run retained **28,690 transcript-cluster features** mapping
to **26,015 distinct Entrez IDs**. Of the 53,617 RMA rows, 24,927 lacked a
single eligible annotation and were excluded with an audit record.

## Outputs and verification

Raw sources are under `data/raw/cit-v1/`. Interim expression is under
`data/interim/cit-v1/`; CSV tables, audits, metadata and validation under
`data/processed/cit-v1/`; and `expression.h5ad` under `data/anndata/cit-v1/`.
H5AD stores samples × transcript clusters in X, clinical data in obs and
annotations in var.

With the complete cache and pinned environment, rebuild without HTTP:

```bash
CIT_RSCRIPT="$PWD/.codex-build/cit-r/bin/Rscript" \
  uv run --offline --locked cit-e-mtab-6128 --offline \
  --output-root /tmp/cit-check/processed \
  --interim-root /tmp/cit-check/interim \
  --anndata-root /tmp/cit-check/anndata
```

[BioStudies deposit](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6128)
· [Original paper](https://www.annalsofoncology.org/article/S0923-7534%2819%2934138-9/fulltext)
