# Frozen cohort input provenance

`taylor-v2-manifest.json` is an unchanged copy of the raw-input manifest used for
the verified Taylor build. It records 58 GEO, cBioPortal and MyGene files/responses,
including original requests, retrieval times, byte sizes and SHA-256 checksums.
It does not contain the expression matrix or clinical/mapping response payloads.

The source snapshot was acquired on November 1, 2026 in America/New_York
(November 2 in the manifest's UTC timestamps). The feature-preserving export uses
these same inputs. The raw expression gzip has SHA-256:

```text
f32c801b0ca2b7df7ea7d9da8c339d9e30fc2c6cf7fbe67e2dae57b1c901b6c7
```

For exact input reproduction, obtain the complete cached `data/raw/taylor-v2/`
directory from a collaborator, then compare its manifest with this reference:

```bash
cmp docs/provenance/taylor-v2-manifest.json data/raw/taylor-v2/manifest.json
uv run --offline --locked taylor-gse21034 --offline --snapshot taylor-v2
```

The pipeline checks each cached file against its manifest before processing it.
Do not put the reference manifest into an empty raw cache: the checksummed
payloads must also be present. A new online acquisition can retrieve changed
upstream results; the manifest records the actual responses, and does not make
remote services immutable. Keep the Git revision and locked dependency file
alongside the raw cache when sharing a reproducible snapshot. Metadata includes
implementation hashes, so code changes can legitimately change output metadata.

## TCGA-PRAD (`tcga-prad-v2`)

`tcga-prad-v2-manifest.json` is an unchanged copy of the raw-input manifest for the
TCGA-PRAD build. It lists 13 files: the cBioPortal datahub RSEM matrix, the clinical
attribute definitions, paged PATIENT/SAMPLE clinical records, the sample inventory, and
the Firehose Legacy patient records used for Gleason. It was acquired on November 11, 2026.
The RSEM matrix (`data_mrna_seq_v2_rsem.txt`, 77,127,361 bytes) is downloaded from cBioPortal datahub
commit `86690e1ed9752b1dcd50b5657f5f05eafa4b6b78` (2026-01-30, the last commit to touch the file), whose
Git LFS pointer records the same hash and size. Its SHA-256:

```text
b7c66db902bb2e7b2e68a1ef8dba4afc6b26f45fb453837bf30163e31d820179
```

This matches the fingerprint recorded by the November 6 pilot run. Every input file is byte-identical
to the earlier `tcga-prad-v1` snapshot, which used the unpinned `master` URL and can no longer be rebuilt
by the current code. To reproduce the
build, get the complete `data/raw/tcga-prad-v2/` directory from a collaborator, then:

```bash
cmp docs/provenance/tcga-prad-v2-manifest.json data/raw/tcga-prad-v2/manifest.json
uv run --offline --locked tcga-prad --offline --snapshot tcga-prad-v2
```

## Cambridge

`cambridge-v2-manifest.json` records the Cambridge GSE70768 Series Matrix,
GPL10558 annotation and the publication supplementary ZIP from Europe PMC
(PMC4588396). It contains requests/hashes, not source payloads. The clinical table
is `mmc2.xlsx` inside the ZIP. See [Cambridge policies](../cambridge.md) for the
subject joins, preservation of GEO values and reconciliation against the paper.
Use the complete `data/raw/cambridge-v2/` directory for offline reproduction.

## Stockholm

`stockholm-v1-manifest.json` records the Stockholm GSE70769 Series Matrix,
GPL10558 annotation and publication clinical supplementary ZIP (PMC4588396).
Source payloads remain under `data/raw/stockholm-v1/`. The publication is shared
with Cambridge, but Stockholm clinical subjects are selected independently by
`STKHLM` identifiers and verified primary-RP specimen type. See the
[Stockholm guide](../stockholm.md) for the two unresolved missing outcomes.

## Mortensen and Long

`mortensen-v1-manifest.json` records GSE46602 expression and GPL570 annotation.
`long-v1-manifest.json` records GSE54460 metadata and the original FPKM supplement.
These snapshots were acquired from empty caches on November 13, 2026; manifests
record original requests, retrieval timestamps, byte sizes and checksums. Payloads
remain under `data/raw/<snapshot>/`. Exact-input offline reproduction requires
those complete raw folders, not just these reference manifests. No raw sequencing
or CEL files are needed by these processed-data adapters.

## CIT

`cit-v1-manifest.json` records the BioStudies inventory, SDRF and all 141
CEL inputs. CEL sizes were verified against the BioStudies inventory, and the
complete payloads were hashed before normalization. The pinned macOS R package
environment is recorded in `env/cit-osx-64.lock`. A raw-cache rebuild also
requires that R environment; see the [CIT guide](../cit.md).
