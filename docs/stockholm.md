# Stockholm GSE70769

The Stockholm adapter uses the same Illumina annotation, missing-expression,
BCR and CSV/H5AD export rules as Cambridge, with explicit Stockholm series,
specimen and subject-selection rules. It downloads the processed Series Matrix,
GPL10558 annotation and the original publication's clinical supplement
(PMC4588396, Supplementary Table 2).

## Run

```bash
uv sync --locked
uv run --locked stockholm-gse70769
```

The default snapshot is `stockholm-v1`. Outputs are separate from Taylor and
Cambridge:

```text
data/raw/stockholm-v1/         Series Matrix, GPL10558, publication ZIP, manifest
data/interim/stockholm-v1/     All deposited probes for selected Stockholm samples
data/processed/stockholm-v1/   expression.csv.gz, genes.csv, clinical.csv,
                              samples.csv, sample_audit.csv,
                              feature_mapping.csv.gz, metadata.json,
                              validation.json, checksums.sha256
data/anndata/stockholm-v1/     expression.h5ad, metadata.json, checksums.sha256
```

H5AD uses samples × probes in `X`, clinical/sample fields in `obs` and gene/RefSeq
annotations in `var`. Original ILMN probe IDs are the unique row/variable keys.
Multiple probes mapped to one gene remain separate. No averaging, imputation,
additional normalization or log transformation is applied.

## Selection and clinical source

Require the GEO primary-RP tumor title and a unique `STKHLM<number>` subject ID.
GEO titles include the spelling `prostatetctomy`; both that deposited form and
`prostatectomy` are recognized explicitly. Stockholm's GEO records do not supply
the Cambridge `sample type` field. The publication must independently confirm
`TUMOUR` and `RP` for each matched subject. Unknown specimen types, duplicate
subjects and unmatched publication identities fail for review.

The publication contains 99 Stockholm subjects, while the deposited expression
matrix contains 94. Only subjects with deposited expression enter this dataset;
five paper-only subjects are not added or assigned expression values. All 94 GEO
profiles remain in the export independently of BCR outcome completeness.

For all selected subjects, use publication BCR status, event time and follow-up
together. Retain original `geo__` characteristics, `paper__` clinical fields,
publication subject IDs and `bcr_clinical_source`. Source event subtypes for the
45 BCR events are `FORMAL (>0.2)`; the adapter does not infer an unreported
confirmation rule or death-handling policy.

Map `Y → 1`, `N → 0`. Events use publication `Time to BCR`; censored observations
use `Total Follow up (months)`. These columns already contain months; do not
convert the separately supplied years columns or apply a second conversion.
Missing outcomes remain nullable, with an explicit eligibility reason. Negative,
nonfinite or malformed times fail. Zero time or an event after recorded follow-up
is retained but flagged ineligible pending source review.

## Verified cohort and unresolved missing outcomes

The November 6, 2026 snapshot contains **94 primary tumors**:

- 45 source-coded BCR events.
- 48 source-coded non-events.
- One missing BCR status.

**92 samples are eligible for BCR analysis: 45 events and 47 censored.** Two
profiles remain in expression, clinical tables and H5AD but are ineligible:

| GEO accession | Subject | Source issue | Eligibility reason |
|---|---|---|---|
| GSM1817928 | STKHLM4753 | Non-event status exists, but total follow-up is UNKNOWN in both GEO and the publication | missing_time |
| GSM1817961 | STKHLM5605 | BCR status and time are unavailable in both sources | missing_status |

These gaps cannot be repaired from the deposited clinical sources. No censoring
time or event code is invented, and these two expression samples are not dropped.

## Probe annotation and missing values

The rules match [Cambridge](cambridge.md): use frozen GPL10558 human Entrez/symbol
annotations, exclude ambiguous or inconsistent gene/symbol identities, then
exclude an entire probe if any selected Stockholm sample has missing expression.
This includes the two samples ineligible for BCR analysis. Every probe's
annotation status, missing count and final selection decision is recorded in
`feature_mapping.csv.gz`.

Of 47,323 deposited probes, seven have 31 missing cells across the 94 samples.
Six otherwise-mapped probes are excluded for missing values. The final dataset
contains **41,615 separate probes representing 29,608 unique Entrez IDs**.
Counts are observations, not hardcoded selection targets. Annotation uses the
deposited platform, not a fresh MyGene/HGNC validation. Upstream processing is
reported as log2 transformation, quantile normalization and ComBat; this adapter
adds no processing of that kind.

## Reproduction and verification

The source manifest records all three input requests, sizes and SHA-256 hashes.
The reference `docs/provenance/stockholm-v1-manifest.json` contains provenance,
not the downloaded payloads. Transfer the entire raw snapshot to reproduce the
same inputs elsewhere, install locked dependencies, then run offline:

```bash
uv run --offline --locked stockholm-gse70769 --offline \
  --output-root /tmp/stockholm-check/processed \
  --interim-root /tmp/stockholm-check/interim \
  --anndata-root /tmp/stockholm-check/anndata
```

Existing differing outputs are protected. Use new output roots after code changes,
or `--snapshot stockholm-v2` for a deliberate upstream-data refresh.

```bash
uv run --locked python -m unittest discover -s tests -v
```

Tests cover both title spellings, explicit cohort separation in the publication,
subject/primary-RP identity checks, missing outcomes, finite-probe selection across
all profiles, preservation of same-gene probes, acquisition, CSV/H5AD round trips,
and byte-identical synthetic offline reproduction. Cambridge and Taylor regression
tests remain included. The shared refactor was also compared against the existing
real Cambridge H5AD: expression, clinical outcomes, sample/probe order and gene
annotations were unchanged. No earlier cohort outputs were overwritten.

Sources: [GSE70769](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70769),
[GPL10558](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558),
[primary publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC4588396/).
