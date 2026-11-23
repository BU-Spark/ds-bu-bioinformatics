# Review resolution and commit scope

The changes add twelve **experimental** cohort workflows alongside eight existing
pipelines, plus a within-cohort PCA/GSVA comparison. Command registration does not
establish full-source validity. No new cohort is advertised as clinically validated.

## Addressed findings

- CancerMap uses `material_type` for tumor/benign status, and source name/organism
  to exclude xenografts. `tissue: Prostate` is an organ field. Donor regions remain audited.
- Supplements require matching known clinical contracts; unresolved origins/censoring
  for deposited times cannot be mixed. Supplemental metadata stays separate from
  unmatched records, and event/duration provenance is explicit.
- Ordinary Cox is blocked for all twelve new experimental cohorts and Sboner's
  outcome-selected extreme groups. Record completeness does not authorize inference.
- Survival screening requires 20 records, 10 events and 5 censorings; this is a minimum
  exploratory safeguard, not a power calculation. Numerical/convergence failures leave
  PCA scores intact. No comparable C-index pairs produce a missing metric and a reason.
- Belfast records published endpoint definitions, with deposited-duration verification
  explicitly false pending full-source reconciliation.
- Exon annotations consider all RefSeq candidates and quarantine source/mapping identity
  conflicts. Unverified endpoint fields are not fabricated across other exon cohorts.
- Sun declares a linear MAS5 scale. Raw count exports retain counts and full-source
  library totals; analysis preprocessing is explicit.
- Failed new local inputs roll back, so corrected source configs/bundles/supplements
  can retry. Shared-builder source-snapshot builds are locked against concurrent updates.
- H5AD `obs`/`var` are independently compared against CSVs; interim source matrices
  are checksummed. Code hashes track dependencies. [Fresh output roots](rebuilding.md)
  are required when old snapshots differ in provenance under current code.

## Deferred rather than silently waived

All twelve new workflows still need complete real-source acquisition, source-schema/count
review and offline-rebuild validation. Public acquisition previously hit HTTP 403 and
DNS failures. DKFZ and CPGEA need legitimately accessible processed files and original
specimen linkage. Exon-array peak memory/disk requirements and Erho's complete supplement
layout remain unmeasured; the dense implementation is **not** claimed operational at full
size. These constraints do not prevent PCA/GSVA on already complete verified exports.

The source-selected Mayo designs and patient overlaps require dedicated analysis methods;
adding complete event times does not remove those restrictions. The PCA/GSVA report
does not contain clinical inference or claim the meeting's eight-cohort totals.

This branch can be committed as experimental pipeline additions and analysis tooling,
with these limitations retained. Full production support for the new cohorts is a
separate validation milestone. Commit and push are authorized by the user; no pull request is created by this workflow.

The shared PCA/GSVA report retains the independent-review Taylor sensitivity caveat:
six module correlations reversed sign under alternate preprocessing. New cohort
workflows remain experimental pending full-source validation.
