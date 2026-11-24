# Rebuilding preserved source snapshots after code changes

Published exports are immutable and compared byte-for-byte. Changes to code provenance,
metadata or H5AD serialization can produce `Existing outputs differ` even when expression
and clinical values are identical. Do not overwrite previous output files or replace
their code hashes with current ones.

Keep the same raw snapshot and use fresh, disjoint output roots:

```bash
uv run cambridge-gse70768 --snapshot cambridge-v2 --offline \
  --output-root data/rebuilt/processed \
  --interim-root data/rebuilt/interim \
  --anndata-root data/rebuilt/anndata
uv run check-outputs --data-dir data/rebuilt --snapshot cambridge-v2
```

The same flags apply to Stockholm, Long and Mortensen. Use each cohort's command and
existing raw snapshot name. This preserves cached biological inputs while recording
the current code. Repeat into a second fresh root to check byte-identical replay under
that code. Original-byte reproduction requires the original code revision, dependencies
and complete source cache. No compatibility shim silently changes provenance.

Code hashes now follow local module dependencies rather than including every unrelated
cohort adapter. This reduces unrelated changes, but real shared-dependency changes must
still change provenance. Independent `check-outputs` verifies source-valued CSV/X,
clinical `obs`, gene `var`, metadata and export checksums. Shared-builder interim matrices
also have their own checksums.

## Running every registered command against real inputs

From a synced Python environment, the validation runner launches the actual cohort
CLIs and independently checks every successful export:

```bash
uv sync --locked
uv run --locked python scripts/run_real_data_validation.py \
  --output data/full-validation --timeout 7200
```

To use the 20 verified public source variants (including the DKFZ processed
subset and experimental CPGEA derivative), first prepare CPGEA as described in
[the RNA-seq guide](authorized_rnaseq.md), then add `--public-sources`. The
validated public snapshots include `cit-fullsource-v1`, `cpgea-public-v1` and
`dkfz-eopc-public-recovery-v1`. These explicit variants preserve their source
limitations and do not claim full controlled raw-deposit coverage.

It records commands, UTC start times, elapsed seconds, return codes and per-cohort
logs in the chosen output root. The default is two smaller cohorts concurrently,
then the large exon-array cohorts sequentially. CIT and the exon-array builders
also share a memory resource lock; run only one large-cohort validation session
at a time. A busy lock is an orchestration conflict, not a source-schema verdict.

Set up CIT's pinned R environment using [its guide](cit.md) first. The runner
detects `.codex-build/cit-r/bin/Rscript` when present. CIT downloads up to four
CELs concurrently by default; `CIT_DOWNLOAD_WORKERS` accepts 1–8. First runs
download several gigabytes and can exceed the per-command two-hour budget on
a slow connection. Validated downloads are retained for retries. DKFZ and CPGEA
require the [documented supplied research inputs](authorized_rnaseq.md); the runner
records their explicit missing-input failures when those have not been supplied.

Rebuild into a second fresh root without HTTP, or choose individual cohorts:

```bash
uv run --locked python scripts/run_real_data_validation.py \
  --output data/full-validation-replay --offline
uv run --locked python scripts/run_real_data_validation.py \
  --output data/full-validation-erho --cohort erho --timeout 7200
```

Inspect both each pipeline return code and its independent checker return code.
Compare file sets and SHA-256 digests between fresh builds. An export passing
these checks establishes source compatibility and reproducibility; normalization,
gene coverage, clinical definitions, sampling and cross-study patient overlap
still determine which analyses are appropriate.
