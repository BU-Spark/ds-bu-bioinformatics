# Real-source processing validation — 2026-11-24

All **20 requested cohorts have processed source-variant outputs**: **4,073 cohort profiles** and **299,190,387 retained feature-by-profile values**. Every snapshot passed its actual source build, CSV/H5AD round trips, checksum checks, and independent output checker. These are expression-processing checks, not clinical validation or approval to pool survival analyses.

The assembled files are under `data/full-validation-all/processed/<snapshot>/` and `data/full-validation-all/anndata/<snapshot>/expression.h5ad`. The former stores features × samples; H5AD stores samples × features with clinical metadata in obs and source/gene annotations in var. The raw snapshots and normalization environment are preserved locally.

[Inventory CSV](results/cohort_inventory_2026-11-24.csv) · [Recorded checks and provenance](provenance/real_source_checks_2026-11-24.json)

| Cohort / source variant | Snapshot | Profiles | Feature rows | Unique gene IDs |
|---|---|---:|---:|---:|
| Belfast / FASTMAN (diagnostic biopsy; radiotherapy + ADT) | `belfast-v1` | 248 | 52,099 | 17,650 |
| Cambridge / CamCap | `cambridge-v2` | 112 | 41,610 | 29,604 |
| CancerMap | `cancermap-v1` | 134 | 12,703 | 12,568 |
| CIT / Kamoun | `cit-fullsource-v1` | 99 | 28,690 | 26,015 |
| CPC-GENE | `cpc-gene-v1` | 213 | 24,414 | 24,414 |
| CPGEA | `cpgea-public-v1` | 134 | 16,034 | 16,031 |
| Den / Thomas Jefferson | `den-v1` | 139 | 259,132 | 15,573 |
| DKFZ / EOPC (public Gerhauser 2018 processed subset) | `dkfz-eopc-public-recovery-v1` | 95 | 16,582 | 16,576 |
| Erho / Mayo | `erho-v1` | 545 | 259,132 | 15,573 |
| Karnes / Mayo | `karnes-v1` | 235 | 259,132 | 15,573 |
| Klein / Cleveland | `klein-v1` | 182 | 4,424 | 286 |
| Long / FFPE RNA-seq | `long-v1` | 100 | 23,280 | 23,206 |
| Mortensen / Aarhus | `mortensen-v1` | 36 | 41,692 | 20,462 |
| Moscow / EMB-RAS | `moscow-v1` | 136 | 36,627 | 36,588 |
| Nakagawa / Mayo (nested case-control; paired DASL panels) | `nakagawa-v1` | 596 | 939 | 886 |
| Sboner / Swedish Watchful Waiting | `sboner-v2` | 281 | 6,144 | 6,100 |
| Stockholm | `stockholm-v1` | 94 | 41,615 | 29,608 |
| Sun / Goodison | `sun-v1` | 79 | 19,695 | 12,427 |
| Taylor/MSKCC | `taylor-v2` | 131 | 25,900 | 19,306 |
| TCGA-PRAD | `tcga-prad-v2` | 484 | 19,088 | 19,064 |

## Source scope and scientific limitations

- **DKFZ:** the public `prostate_dkfz_2018` Gerhauser 2018 matrix contains 118 tumor regions from 95 early-onset ICGC donors. Exact sample/patient links are verified; one region per donor is retained by original sample ID and the other 23 are audited. The EGA description lists 96 donors plus nine controls: the public subset is not the complete controlled raw deposit. The cBioPortal profile description has incorrect neuroblastoma/sample totals; this discrepancy is recorded. Original gene symbols are mapped through cached exact human responses. The follow-up column lacks an explicit unit, so no survival durations are assigned.
- **CPGEA:** an experimental public secondary derivative provides 134 profiles. The source clinical IDs match, but ten FOXA1 group labels disagree with the clinical supplement; both are preserved. Unsupported XLOC rows remain in the preparation audit. A source `Dead` label is not converted to recurrence or nonrecurrence. Duration semantics remain unresolved.
- **CIT:** 141 original CELs were acquired. The pinned standard reader rejects OCT_069 and OCT_056; an alternative reader returns extreme raw intensities for OCT_069. Both source files are preserved and their two tumor-array exclusions are audited. RMA used 139 readable arrays, including 40 nontumor arrays; the final export has 99 tumors and 75 BCR-eligible records with 13 events.
- **Klein:** the deposited processed matrix has only 18,846 PSR rows, primarily chromosome 1; the export retains 4,424 rows representing 286 genes. Whole-transcriptome recovery requires a separately validated raw-CEL workflow.
- **Nakagawa:** the conservative source/RefSeq policy retains 939 rows representing 886 genes. Targeted panel coverage and unmapped GenBank-only targets limit broad module analyses.
- **Endpoints and sampling:** many cohorts have status without verified individual duration; those profiles remain in expression and are ineligible for survival fits. The experimental analysis restrictions remain in metadata. Sboner outcome-selected groups, Mayo/Decipher overlap, treatment differences and case-cohort sampling still require appropriate analysis designs. Eligibility counts express record completeness, not design approval.
- **Cross-study identity and scale:** totals are cohort profiles, not globally deduplicated people. Expression scales differ (raw counts, RSEM/log2, FPKM/RPKM/TPM, microarrays); no automatic pooling or cross-study normalization was performed. Gene-level and probe/PSR-level rows remain distinct.

## Reproducibility and operation

All 159 tests pass. Thirteen full real-source cohorts additionally reproduced byte for byte in two fresh offline builds. The other real source variants passed structural and value round trips; full second-build byte comparisons for all of those large snapshots have not been claimed. Cached annotation reuse requires identical URL, method, request body and verified response bytes, and preserves original acquisition metadata.

See [public/authorized RNA-seq routes](authorized_rnaseq.md), [CIT setup and QC](cit.md), and [running/rebuilding](rebuilding.md). Prepare the public CPGEA inputs and pinned CIT environment, then use `python scripts/run_real_data_validation.py --public-sources` for a fresh output root. Rebuild requests after code changes require fresh output directories; changing biological input requests requires a new raw snapshot.

The raw and processed matrices are intentionally excluded from Git. The versioned inventory and checker/provenance record contain the observed counts and validation evidence.
