# Sun / Goodison — GSE25136

Real-source acquisition and independent CSV/H5AD checks pass for **79 profiles**.
This cohort is **not survival-ready** because individual follow-up durations have
not been verified. See the [full-source report](full_source_validation.md).

```bash
uv run --locked sun-gse25136 --snapshot sun-v1
uv run --locked check-outputs --snapshot sun-v1
```

Module execution without refreshing installed entry points:
`PYTHONPATH=src .venv/bin/python -m ds_bu_bioinformatics.sun_pipeline`.

## Inputs and policies

- [GSE25136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25136): GEO describes
  79 primary tumors, including 39 recurrence cases. These are reference counts;
  this implementation has not yet verified them from a complete downloaded matrix.
- [GSM617640](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM617640): example
  `Prostate cancer primary tumor PR58`, `recurrence status: Non-Recurrent`;
  processing is MAS5 with global scaling.
- [GPL96](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96): frozen platform
  Entrez, symbol, species and RefSeq annotation.

The adapter requires source primary-tumor tissue and unique PR* title-derived subject
labels; it retains missing outcomes. Duplicate patients/assays, unknown tissue or
platform/schema changes fail. No claim of independent registry linkage is made.

Original numeric expression values and separate probe rows are preserved. Only
single human Entrez/symbol mappings are selected, excluding ambiguous/conflicting
annotations. Features missing in any selected subject are excluded and audited;
values are never imputed or averaged. MAS5 signal is not automatically relabeled
as log2 expression. The existing module scorer's maximum-value log heuristic remains
provisional; normalization needs review before biological interpretation.

## Recurrence status is not survival time

`Recurrent` / `Non-Recurrent` maps to `recurrence_event` 1 / 0. Missing status stays
missing. The class is **not relabeled BCR, DFS or OS**. Every `recurrence_months` is
null, `recurrence_eligible` is false, and eligibility reasons explain missing duration
or status. The common `_time_unit = months` indicates the export field's intended unit;
it does not assert that a duration exists.

PCA scores are still produced. Cox fitting, HR and C-index are skipped with a note;
a run containing only cohorts without usable survival data still writes a report
and a plot explaining why no survival models were fitted. A structural output-check
PASS is not evidence of survival readiness.

To enable survival analysis, obtain a patient-linked event/censoring-time table and
verify the PR-label mapping, outcome definition, time origin, units and ascertainment.
Never assign median follow-up, zero, or a guessed duration to individual patients.
The [clinical supplement interface](clinical_linkage.md) can import verified times
through exact exported sample IDs, reject conflicting statuses and replay offline.
Use a new snapshot for supplemented data. This changes eligibility only for matched
records and preserves the original recurrence labels and missing source dates.

## First full-source verification

Download and inspect all source schemas, subject labels and recurrence codes; compare
counts and feature exclusions with the paper. Run `check-outputs`, then rebuild into
separate output roots with `--offline` and compare all file hashes, as described in
[the Belfast guide](belfast.md). No full-cohort reference manifest is supplied until
real acquisition and reproducibility have actually been verified.
