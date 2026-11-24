# Remaining source-data validation and clinical work

All 20 requested cohorts now have processed public source variants with verified
CSV/H5AD exports (4,073 cohort profiles). The twelve new workflows remain
experimental for clinical interpretation and analysis design. See the
[full-source validation report](full_source_validation.md) and
[coverage table](pipeline_coverage.md). This file records remaining scientific
requirements rather than unimplemented acquisition commands.

DKFZ has a verified public Gerhauser 2018 processed subset (95 donors); the complete
controlled raw deposit is not acquired. CPGEA uses an explicitly documented public
secondary derivative (134 profiles). Both still support supplied processed sources.

| Cohort | Implemented workflow | Remaining source/analysis issue and next step |
|---|---|---|
| CancerMap | Native GSE94767/GPL5175 acquisition; malignant human donor/region selection; source-feature mapping | Full deposited-source run passed (134 donors); review region sensitivity. Locate and verify linked BCR clinical data. |
| CPC-GENE | Native processed TSV plus both metadata platforms; exact CPCG/GSM join; source batch/replicate policies | Complete 213-profile matrix passed; verify overlap with clinical supplement IDs. Source probe summarization, batch correction and CPCG0398 replicate averaging are preserved, not undone. |
| DKFZ / EOPC | Dedicated authorized processed-RNA workflow with donor/area mapping | Public Gerhauser 2018 subset passed for 95 donors; full controlled-deposit coverage and duration-unit reconciliation remain separate. |
| Erho / Mayo | Native processed PSR/GEO workflow and metastatic-event labels | Full matrix can be large. Obtain/adjudicate individual event/censoring times; assess selected sampling and overlap. |
| Nakagawa / Mayo | Native exact paired-panel workflow, source clinical fields, RefSeq mapping and PCSM class | Full paired-panel export passed for 596 profiles. Verify PCSM duration separately; do not infer final systemic-event status from SYS/PSA/NED groups. GenBank-only/Affymetrix target mappings currently remain excluded/audited. |
| Karnes / Mayo | Native GSE62116/GPL5188 workflow | Verify individual clinical event/time linkage, case-cohort sampling and Mayo patient overlap. |
| Den / Thomas Jefferson | Native GSE72291/GPL5188 workflow with RT context | Obtain linked BCR/metastasis durations and appropriate postoperative-radiotherapy time origins. |
| Klein / Cleveland | Native GSE62667/GPL5188 workflow | Obtain verified event/time linkage; preserve the source's SCAN processing and study context. |
| CPGEA | Dedicated processed-RNA source workflow | Public 134-profile derivative passed; resolve source-group discrepancies and clinical duration semantics. |
| Moscow / EMB-RAS | Native submitted mRNA counts, exact assay joins, tumor/mRNA selection and BCR class | Complete 136-profile source export passed; obtain verified individual follow-up. Raw export counts remain unchanged; library totals support separately reviewed downstream normalization. |
| Sun / Goodison | Native GSE25136/GPL96 workflow and recurrence class | Need verified individual recurrence/censoring times for survival. MAS5 scale/normalization requires analysis review. |
| Belfast / FASTMAN | Native GSE116918/GPL25318 workflow, separate BCR/metastasis durations | Full 248-profile export passed; reconcile event counts and endpoint time origins against clinical methods; retain biopsy/RT/ADT context. |

## Available mechanisms to resolve gaps

- [Clinical supplement import](clinical_linkage.md) fills verified missing values using
  exact exported sample IDs, rejects contradictions and caches all input provenance.
- [DKFZ/CPGEA source configuration](authorized_rnaseq.md) maps original source columns,
  selects primary tumors, handles multiple areas explicitly and generates all exports.
- [Generic curated import](cohort_import.md) is an optional route for previously prepared
  research inputs, with full source/audit/hash requirements.
- All new workflows retain nullable outcomes and support expression/PCA even when
  survival fitting is unavailable. None infers individual dates from group labels,
  literature totals or median follow-up.

## Source leads for unresolved clinical links

- CancerMap: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94767),
  [DESNT primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5669460/).
- CPC-GENE: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107299),
  [clinical Table S1](https://pmc.ncbi.nlm.nih.gov/articles/PMC6511374/),
  and [Fraser et al. clinical Supplementary Data 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC8556363/).
  Their overlap and sample crosswalks with the expression matrix still need verification.
- Nakagawa: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10645),
  [example clinical record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM268556),
  [primary study design](https://pmc.ncbi.nlm.nih.gov/articles/PMC2565588/).
- Erho: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46691).
- Karnes: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62116).
- Den: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72291).
- Klein: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62667).
- Moscow: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229904).
- DKFZ: [EGA source/access record](https://ega-archive.org/datasets/EGAD00001004791).
- CPGEA: [primary publication](https://www.nature.com/articles/s41586-020-2135-x),
  [GSA-Human record](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000099).
