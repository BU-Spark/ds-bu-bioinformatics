# Frozen Taylor input provenance

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

## TCGA-PRAD (`tcga-prad-v1`)

`tcga-prad-v1-manifest.json` is an unchanged copy of the raw-input manifest for the
TCGA-PRAD build. It lists 13 files: the cBioPortal datahub RSEM matrix, the clinical
attribute definitions, paged PATIENT/SAMPLE clinical records, the sample inventory, and
the Firehose Legacy patient records used for Gleason. It was acquired on November 11, 2026.
The RSEM matrix (`data_mrna_seq_v2_rsem.txt`, 77,127,361 bytes) has SHA-256:

```text
b7c66db902bb2e7b2e68a1ef8dba4afc6b26f45fb453837bf30163e31d820179
```

This matches the fingerprint recorded by the November 6 pilot run. To reproduce the
build, get the complete `data/raw/tcga-prad-v1/` directory from a collaborator, then:

```bash
cmp docs/provenance/tcga-prad-v1-manifest.json data/raw/tcga-prad-v1/manifest.json
uv run --offline --locked tcga-prad --offline --snapshot tcga-prad-v1
```
