# Adding verified clinical follow-up without changing deposited source fields

The shared cohort builder and `import-cohort` accept an explicit clinical supplement.
This provides a reproducible route for missing clinical links; it **does not locate,
verify or invent unavailable dates**. It applies to Belfast, Sun, Nakagawa,
Mortensen, Long, CancerMap, CPC-GENE, Erho, Karnes, Den, Klein, Moscow, the
authorized DKFZ/EOPC and CPGEA workflows, and imported bundles. Taylor, TCGA,
Cambridge, Stockholm, CIT and Sboner retain their existing separate builders.

## Contract

Prepare a CSV with one row per **selected output sample and declared endpoint**:

```csv
sample_id,endpoint,event,time,time_unit
SAMPLE_ID_FROM_EXPORTED_SAMPLES_CSV,recurrence,1,INDIVIDUAL_EVENT_TIME,months
ANOTHER_EXPORTED_SAMPLE_ID,recurrence,0,INDIVIDUAL_CENSORING_TIME,years
```

These are placeholders, not clinical data. Use exact `samples.csv` sample IDs.
Nakagawa uses `GSE10645:patient:<source numeric label>` for a paired sample. Preserve
leading zeroes. If original clinical records use another identifier, establish and
review that crosswalk before supplying the CSV; the importer does not fuzzy-match
patient names, drop unmatched records or infer aliases.

- `event`: exactly `0` or `1`, with the endpoint's actual meaning.
- `time`: a finite nonnegative individual duration. Zero remains ineligible.
- `time_unit`: `months` or `years`; years convert to months by multiplying by 12.
- `endpoint`: an existing prefix from the cohort. Sun uses `recurrence`, Belfast
  `bcr`/`metastasis`, Nakagawa `pcsm`. Endpoint relabeling requires a separate reviewed adapter change.

A JSON provenance file must identify `study_accession`, the clinical source, and each
endpoint's `definition`, `time_origin`, `censoring` and `sample_linkage`. Start from
[the provenance template](templates/clinical_linkage_provenance.example.json) and
replace its placeholders with the actual evidence. Describing a clinical source is
not a substitute for verifying it.

## Run and replay

```bash
uv run --locked sun-gse25136 --snapshot sun-linked-v1 \
  --clinical-file /path/to/verified_clinical.csv \
  --clinical-provenance /path/to/verified_provenance.json
uv run --locked check-outputs --snapshot sun-linked-v1
```

Use a new snapshot when adding/changing a supplement. The table and provenance are
copied into `data/raw/<snapshot>/clinical-linkage/` and frozen in the input manifest.
Outputs contain the full submitted table in `clinical_linkage_audit.csv`, linked
values in separate `linked__*` columns, original source tokens unchanged, and effective
endpoint provenance in metadata. Unmatched selected samples retain their source outcomes.

For an offline rebuild, use the same snapshot, omit the original clinical file paths,
and supply separate output roots. The cached supplement is replayed without network
access or requiring the original files:

```bash
uv run --locked sun-gse25136 --snapshot sun-linked-v1 --offline \
  --output-root data/rebuilt/processed --interim-root data/rebuilt/interim \
  --anndata-root data/rebuilt/anndata
```

Unknown sample IDs, duplicate sample/endpoint pairs, undeclared endpoints, missing
provenance, malformed dates and conflicts with known deposited event/time values fail.
This tool fills verified missing data; it does not silently resolve contradictory
clinical sources. For contradictions, obtain/adjudicate an authoritative source and
implement an explicitly documented cohort-specific policy.

## Definition compatibility and analysis permission

Supply the exact documented source `definition`, `time_origin` and `censoring` strings
where known. Numeric agreement alone is insufficient. If source times exist but their
origin or censoring is unresolved, the supplement is rejected until a reviewed source
contract is encoded. The original endpoint metadata stays unchanged; supplemental
definitions apply only to matched rows, with separate event and duration provenance.
Rejected new local inputs are rolled back so corrected inputs can be retried in the
same source snapshot. Accepted inputs remain immutable.

`*_eligible` marks complete positive-time records, not approval for ordinary Cox.
All twelve experimental cohorts and Sboner are blocked from ordinary Cox in module-scores.
