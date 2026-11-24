# Nakagawa / Mayo — GSE10645

The complete paired-panel source run passes independent CSV/H5AD checks for
**596 profiles**, retaining **939 rows / 886 mapped genes** under the conservative
RefSeq policy. Targeted coverage, unmapped GenBank-only targets and clinical
duration interpretation remain limitations. See the [full-source report](full_source_validation.md).

```bash
uv run --locked nakagawa-gse10645 --snapshot nakagawa-v1
uv run --locked check-outputs --snapshot nakagawa-v1
```

The source series lists 1,192 assay records across GPL5858 and GPL5873 for 596 patients.
The adapter reads both deposited Series Matrices and platforms, requires exactly the
same patient labels on both panels and identical source clinical covariates, then
joins the assays by the title-derived numeric patient label. No positional pairing,
subject dropping or replicate averaging is allowed. A selected sample is a **composite
of two panel assays** with ID `GSE10645:patient:<label>`; both source GSM IDs are retained.
`assay_linkage.csv` records every contributing original assay.

Feature IDs are `<platform>:<original target ID>`. The deposited platforms are
**gene-target representations that already average probes**. The adapter adds no
aggregation, but cannot recover original raw probes from these matrices. Original
numeric signals remain unchanged in CSV/H5AD.

Both GPL tables expose `ID`, `SequenceSource`, and `GB_ACC`. RefSeq targets are mapped
through the existing checksummed MyGene cache with human identity and returned-RefSeq
validation. GenBank-only/Affymetrix-based targets are currently excluded and audited,
not guessed from their labels. Mapping candidates are saved separately in
`gene_mapping_candidates.csv.gz`. Broadening annotation is possible after an explicit
validated accession-mapping extension; limited targeted-gene/module coverage remains.

## Clinical fields are not interchangeable

The source includes PCa-specific death labels, years of follow-up after RRP, PSA-rise
times from RRP and an interval from PSA recurrence to systemic disease. These are
retained; years convert to months in separate `source_*` fields. The adapter does not
add unlike time origins together or interpret a missing systemic interval as a proven
non-event. The primary study documents later systemic progression in some PSA controls,
so **SYS/PSA/NED group membership is not the final systemic-event indicator**.

`pcsm_event` preserves the deposited Yes/No death label. `pcsm_months` remains missing
until an authoritative patient-linked event/censoring duration is supplied. The
`Years of Followup after RRP` label alone is not silently treated as a separately
verified time of death. PCA can run; Cox fitting is skipped when no eligible PCSM events
exist. This also avoids interpreting the study's nested case-control sample as an
unselected population cohort.

To enable PCSM survival analysis, verify that an authoritative duration corresponds
to each patient's PCa death or appropriate censoring time, document the time origin
and competing-event policy, and use [the clinical supplement interface](clinical_linkage.md)
in a new snapshot. If source-panel clinical fields disagree, adjudicate the discrepancy;
the pipeline deliberately fails instead of choosing one panel's outcome.

## Sources

- [GSE10645](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10645)
- [GPL5858](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL5858)
- [GPL5873](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL5873)
- [Example clinical fields, GSM268556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM268556)
- [Primary paper and study design](https://pmc.ncbi.nlm.nih.gov/articles/PMC2565588/)

On the first complete-source run, inspect all subject labels, panel clinical agreement,
annotation exclusions and group/death counts, then rebuild offline into separate roots
and compare all file hashes. Do not relax source validation to force literature counts.
