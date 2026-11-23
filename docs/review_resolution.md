# Dataset pipeline scope and validation

This change adds twelve **experimental** cohort workflows alongside eight existing
pipelines. Command registration does not establish full-source validity. No new
cohort is advertised as clinically validated. Module-scoring algorithms, gene-module
tables, PCA/GSVA tooling and analysis results are outside this PR.

## Dataset correctness and reproducibility

- CancerMap uses `material_type` for tumor/benign status and source name/organism
  to exclude xenografts. `tissue: Prostate` is an organ field. Donor regions remain audited.
- Supplements require matching known clinical contracts. Unresolved origins/censoring
  for deposited times cannot be mixed. Supplemental metadata stays separate from
  unmatched records, with explicit event/duration provenance.
- Exports distinguish complete records from analysis permission through restrictive
  `analysis_policy` metadata. Downstream consumers must enforce those restrictions.
- Belfast records published endpoint definitions, with deposited-duration verification
  explicitly false pending full-source reconciliation.
- Exon annotations consider all RefSeq candidates and quarantine source/mapping identity
  conflicts. Unverified endpoint fields are not fabricated across other exon cohorts.
- Sun declares a linear MAS5 scale. Raw-count exports retain counts and full-source
  library totals for separately reviewed downstream normalization.
- Failed new local inputs roll back, allowing corrected configs/bundles/supplements
  to retry. Shared-builder source snapshots are locked against concurrent updates.
- H5AD `obs`/`var` are independently compared against CSVs; interim source matrices
  are checksummed. Code hashes track dependencies. [Fresh output roots](rebuilding.md)
  are required when existing snapshots differ in provenance under current code.

## Remaining validation

All twelve new workflows still need complete real-source acquisition, source-schema/count
review and offline-rebuild validation. Public acquisition previously hit HTTP 403 and
DNS failures. DKFZ and CPGEA need legitimately accessible processed files and original
specimen linkage. Exon-array peak memory/disk requirements and Erho's complete supplement
layout remain unmeasured; the dense implementation is not claimed operational at full size.

Selected Mayo study designs and patient overlap require dedicated analysis review;
adding complete event times does not remove those restrictions. Full production support
for the new cohorts remains a separate validation milestone.
