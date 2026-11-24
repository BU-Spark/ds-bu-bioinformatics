# Dedicated DKFZ/EOPC and CPGEA RNA-seq pipelines

Both studies now have dataset-specific source workflows:

| Cohort | Command | Required accession | Default snapshot |
|---|---|---|---|
| DKFZ / EOPC | `dkfz-eopc` | EGAD00001004791 | `dkfz-eopc-v1` |
| CPGEA | `cpgea` | PRJCA001124 | `cpgea-v1` |

These workflows accept processed expression plus the original sample/clinical source
table and a documented column map. They perform identity checks, tumor selection,
patient/region handling, gene annotation, endpoint construction, library normalization
metadata, checksummed caching and CSV/AnnData export. They do **not** require a prebuilt
generic cohort bundle or manually prepared gene/sample/export files.

They require legitimately accessible research files. They do not acquire controlled
FASTQs, implement read alignment, grant EGA/GSA access, or claim that such files have
been downloaded here. Synthetic source-format integration and offline-rebuild tests
cover both workflows; real-source files and clinical semantics remain unverified.

## Input files

1. A processed gene-by-assay **TSV or TSV.gz**, with original gene identifiers in its
   first column and unique source assay IDs in all other columns. Keep original numeric
   values. Supported consistent identifier namespaces: Entrez IDs, Ensembl gene IDs
   (including version suffixes), or gene symbols; exact human mappings are cached.
2. The original sample/clinical **TSV or TSV.gz** containing actual assay, patient,
   specimen and available endpoint columns. Preserve string IDs and leading zeroes.
3. A JSON source configuration, with the exact accession, source citation, expression
   kind/scale, source column names, specimen value definitions and endpoint semantics.

Examples are provided for [DKFZ](templates/dkfz_source_config.example.json) and
[CPGEA](templates/cpgea_source_config.example.json). They contain placeholder column
names, specimen labels and citations; edit them to match the actual authorized sources.
They are **not assertions about the original study's exact file schema**.

`expression_kind` is `raw_counts`, `tpm` or `fpkm`. Original values stay unchanged in
exports. Raw counts must be complete, nonnegative integers with positive library totals;
library totals are frozen before annotation filtering. Normalization of raw counts,
TPM or FPKM for downstream analysis is outside this pipeline change.

`sample_columns` maps output concepts to actual source table columns. Every expression
column must map to exactly one metadata `expression_id`; no partial-string or positional
matching is performed. All specimen values must be explicitly listed as tumors or
excluded tissue. Source metadata profiles absent from expression are retained in the
audit, with that exclusion reason.

The default `replicate_policy` is `error` when multiple primary areas/assays belong to
one patient. A deliberately chosen `lowest_sample_id` policy retains one original region
independently of outcomes, records the others, and never averages. Conflicting known
clinical outcomes/times across areas fail even with a selection policy.

Each endpoint configuration declares:

- `definition`, `time_origin` and `censoring`;
- actual `event_column` and explicit `event_codes` mapping to 0/1;
- actual `time_column` and `time_unit` (`months` or `years`).

If a field is unavailable, set its column to JSON `null`. Missing outcomes/times stay
missing and ineligible. A missing date does not block expression or PCA. No source
study total, group label or median follow-up is substituted for an individual outcome.

## First run

```bash
uv run --locked dkfz-eopc \
  --expression-file /path/to/authorized_expression.tsv.gz \
  --sample-file /path/to/original_clinical.tsv \
  --source-config /path/to/verified_dkfz_source_config.json

uv run --locked cpgea \
  --expression-file /path/to/processed_expression.tsv.gz \
  --sample-file /path/to/original_clinical.tsv \
  --source-config /path/to/verified_cpgea_source_config.json
```

The first run freezes all three files and mapping responses in the raw snapshot.
For a cached offline rebuild, the original files need not remain available:

```bash
uv run --locked dkfz-eopc --offline \
  --output-root data/rebuilt/processed --interim-root data/rebuilt/interim \
  --anndata-root data/rebuilt/anndata
```

Then use `check-outputs`, inspect sample/feature/endpoint audits, compare rebuild hashes,
and separately establish clinical suitability. Changing inputs or their interpretation
requires a new snapshot. Both commands also accept the optional
[verified clinical supplement](clinical_linkage.md) interface.

## Source-specific constraints

DKFZ must use the EOPC accession listed above; a different convenient DKFZ/cBioPortal
study is not interchangeable. Multi-area samples need explicit donor/area mapping.
See [the EGA source/access record](https://ega-archive.org/datasets/EGAD00001004791).

### Public processed DKFZ alternative

The correct Gerhauser 2018 study also has a public processed cBioPortal deposit,
`prostate_dkfz_2018` (PMID 30537516). It is distinct from `prad_guru_2017`.
It can be processed without obtaining controlled raw EGA reads:

```bash
uv run --locked dkfz-eopc --public-source --snapshot dkfz-eopc-public-recovery-v1
```

This freezes cBioPortal/datahub commit
`86690e1ed9752b1dcd50b5657f5f05eafa4b6b78`; the verified snapshot is
`dkfz-eopc-public-recovery-v1` (`dkfz-eopc-public-v1` is the command's default
for a fresh first run).
The actual expression matrix contains 118 primary-tumor regions from 95 donors,
all linked exactly to ICGC / EO / Tumor / localized source clinical records.
One region per donor is selected by the lowest original sample ID, independently
of outcomes; the other 23 regions are audited. This is the available public
processed subset, not a claim to have acquired all 96 donors and nine controls
in the controlled EGA description.

The source Entrez column is empty. Original Hugo symbols and row-qualified
duplicate feature keys are preserved; exact human symbol-to-Entrez mappings are
cached, and ambiguous/retired labels are excluded with an audit. RPKM-labelled
source values are retained unchanged. The profile's descriptive text incorrectly
mentions neuroblastoma and mismatched totals; the adapter verifies the study ID,
publication and actual RNA/sample/clinical identities and retains that limitation.
Deposited BCR status is retained. The numeric follow-up column has no explicit
unit, so it remains verbatim source data and is not converted into eligible
survival durations.

For CPGEA, do not assume the total 208-patient study size equals the expression cohort:
the publication reports RNA-seq for a subset. Use the actual expression columns and
verified clinical crosswalk. See [the primary publication](https://www.nature.com/articles/s41586-020-2135-x)
and [GSA-Human record](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000099).

### Experimental public CPGEA derivative

A secondary public repository also provides the 134 expression profiles and the
original study clinical supplement. Prepare its source inputs, then run the
existing CPGEA pipeline:

```bash
uv run --locked python scripts/prepare_public_cpgea.py
uv run --locked cpgea --snapshot cpgea-public-v1 \
  --expression-file data/raw/cpgea-public-discovery/prepared/expression.tsv \
  --sample-file data/raw/cpgea-public-discovery/prepared/samples.tsv \
  --source-config data/raw/cpgea-public-discovery/prepared/source-config.json
```

The helper freezes `nationstrong/CPGEA` commit
`13a8137dce2e2f447702354502d4b903dec40080`, checks exact source identities,
combines the two deposited FOXA1 analysis groups without averaging, and preserves
TPM values. The filenames' `normal` group is a FOXA1 comparison group, not a
tumor/normal tissue pairing. All 134 sample labels match the original clinical
sequencing IDs. Ten group labels disagree with the clinical FOXA1 field; both
are retained without adjudication. All 1,998 unsupported XLOC/non-Ensembl rows
remain in the original files and preparation audit. `Dead` in the recurrence
column remains an unknown recurrence event. Duration semantics remain unresolved,
and no survival durations are assigned. This derivative remains experimental.

After preparing CPGEA and setting up CIT's R environment, the real-data runner
can choose both explicit public source variants:

```bash
uv run --locked python scripts/run_real_data_validation.py --public-sources \
  --output data/full-validation-public
```
