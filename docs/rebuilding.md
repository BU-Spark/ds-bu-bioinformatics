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
