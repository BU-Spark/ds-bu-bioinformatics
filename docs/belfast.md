# Belfast / FASTMAN — GSE116918

New adapter: **implemented and tested with synthetic source-shaped fixtures; full-cohort
acquisition, counts and reproducibility have not yet been verified**. The development
execution environment could browse indexed GEO records but could not resolve GEO's
bulk-download host. Do not label this cohort analysis-ready until the checks below pass.

## Run

```bash
uv run --locked belfast-gse116918 --snapshot belfast-v1
uv run --locked check-outputs --snapshot belfast-v1
uv run --locked module-scores --snapshot belfast-v1
```

Without installing the refreshed CLI entry points, the equivalent pipeline command is
`PYTHONPATH=src .venv/bin/python -m ds_bu_bioinformatics.belfast_pipeline`.

The pipeline uses the shared checksummed input cache and atomic CSV/AnnData export.
It fetches the GSE116918 Series Matrix and full GPL25318 SOFT annotation; it does not
need raw CEL files. Distinct probes remain distinct. No additional expression
normalization, imputation, gene averaging or clinical time inference is performed.

## Source and selection

- [GSE116918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116918): diagnostic
  prostate biopsies from a cohort treated with radical radiotherapy and ADT. GEO lists
  248 samples; this is a source reference, **not a verified export count**.
- [Example sample GSM3263779](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3263779)
  provides source tissue, S1236-* title, BCR/metastasis status and separate durations.
- [GPL25318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL25318) exposes
  `ID`, `Ensembl Gene ID`, `Gene Symbol`, and `Entrez Gene` columns.

Require primary-prostate-cancer tissue and unique S1236-* title-derived subject labels.
Unexpected specimen types, platform changes or duplicate subjects fail for review;
there is no arbitrary replicate selection. Source labels are not an independently
validated patient registry. Missing clinical data do not remove expression samples.

Use single Ensembl/symbol annotations, retaining original Entrez labels separately.
Ambiguous, unmapped and conflicting annotations are quarantined in the feature audit.
Deposited gene labels can be outdated; no external annotation service is queried.

## Endpoints

| Export prefix | Event source field | Duration source field |
|---|---|---|
| `bcr` (primary) | `bcr event (1=yes, 0=no)` | `follow-up time (bcr, months)` |
| `metastasis` | `met event (1=yes, 0=no)` | `follow-up time (met, months)` |

Each exports `_event`, `_months`, `_eligible`, `_eligibility_reason`, and `_time_unit`.
Use each endpoint's own deposited duration; never substitute metastasis follow-up
for BCR or vice versa. Missing/zero times remain nullable/zero and ineligible. Negative,
nonfinite, malformed times or unknown status codes fail. GEO `UNKNOWN` is treated
as missing. Treatment context is cohort-level, not a fabricated patient-specific label.

`module-scores` currently tests the primary endpoint, BCR. Exporting metastasis does
not imply a second metastasis Cox model is fitted automatically.

## Required first real-data checks

1. Download and freeze sources; inspect any schema/selection failure instead of relaxing
   validation to reach an expected cohort count.
2. Inspect `validation.json`, sample/feature audits and the eligible/source-coded counts
   for **both** endpoints. Reconcile any differences with the study's published counts.
3. Confirm endpoint time origin/meaning against the primary clinical methods before
   interpreting survival associations. Retain RT/ADT context in comparisons.
4. Rebuild offline into separate roots and compare every output hash:

```bash
uv run --locked belfast-gse116918 --offline --snapshot belfast-v1 \
  --output-root data/rebuilt/processed --interim-root data/rebuilt/interim \
  --anndata-root data/rebuilt/anndata
```

Do not compare the frozen snapshot to outputs generated after code changes without
reviewing provenance: code hashes are intentionally part of metadata.
