# Pipeline coverage for all 20 requested cohorts

The working tree registers a dataset-specific command for every listed cohort.
The twelve new workflows remain **experimental** for clinical interpretation and
analysis design. All 20 requested cohorts now have real-source processed outputs,
using explicit public processed variants for DKFZ and CPGEA. Those variants do not
claim equivalence to controlled raw-sequencing deposits or resolve missing clinical
semantics. The generic `import-cohort` utility is optional and is
**not** counted as a dataset-specific pipeline.

Real-source builds and independent CSV/H5AD/checksum checks now cover **20 source
variants and 4,073 cohort profiles**. All 159 tests pass; 13 full real cohorts also
reproduced byte for byte offline. Initial network failures were overcome through
completed acquisitions and exact verified cache reuse. See the
[full-source validation report](full_source_validation.md) for actual counts,
source variants, QC exclusions and remaining scientific limitations.

| Cohort | Accession | Dataset-specific command | Acquisition / guide |
|---|---|---|---|
| TCGA-PRAD | TCGA-PRAD | `tcga-prad` | Existing cBioPortal workflow; [dictionary](tcga_prad_data_dictionary.md) |
| Taylor / MSKCC | GSE21034 | `taylor-gse21034` | Existing GEO + cBioPortal workflow; [dictionary](data_dictionary.md) |
| Cambridge / CamCap | GSE70768 | `cambridge-gse70768` | Existing public workflow; [guide](cambridge.md) |
| Stockholm | GSE70769 | `stockholm-gse70769` | Existing public workflow; [guide](stockholm.md) |
| CancerMap | GSE94767 | `cancermap-gse94767` | Native public GEO; [guide](additional_geo.md) |
| CIT / Kamoun | E-MTAB-6128 | `cit-e-mtab-6128` | Existing BioStudies/CEL-RMA workflow; [guide](cit.md) |
| CPC-GENE | GSE107299 | `cpc-gene-gse107299` | Native processed matrix + two-platform metadata; [guide](additional_geo.md) |
| Long | GSE54460 | `long-gse54460` | Existing GEO/FPKM workflow; [guide](long.md) |
| DKFZ / EOPC | EGAD00001004791 | `dkfz-eopc` | Dedicated authorized processed-source workflow; [guide](authorized_rnaseq.md) |
| Belfast / FASTMAN | GSE116918 | `belfast-gse116918` | Native public GEO; [guide](belfast.md) |
| Sboner / Swedish WW | GSE16560 | `sboner-gse16560` | Existing public workflow; see README |
| Erho / Mayo | GSE46691 | `erho-gse46691` | Native public PSR matrix/GEO; [guide](additional_geo.md) |
| Nakagawa / Mayo | GSE10645 | `nakagawa-gse10645` | Native paired DASL panels; [guide](nakagawa.md) |
| Karnes / Mayo | GSE62116 | `karnes-gse62116` | Native public GEO; [guide](additional_geo.md) |
| Den / Thomas Jefferson | GSE72291 | `den-gse72291` | Native public GEO; [guide](additional_geo.md) |
| Klein / Cleveland | GSE62667 | `klein-gse62667` | Native public GEO; [guide](additional_geo.md) |
| CPGEA | PRJCA001124 | `cpgea` | Dedicated processed-source workflow; [guide](authorized_rnaseq.md) |
| Moscow / EMB-RAS | GSE229904 | `moscow-gse229904` | Native public mRNA count matrix/GEO; [guide](additional_geo.md) |
| Mortensen / Aarhus | GSE46602 | `mortensen-gse46602` | Existing public workflow; [guide](mortensen.md) |
| Sun / Goodison | GSE25136 | `sun-gse25136` | Native public GEO; [guide](sun.md) |

## What is still required

Command registration is complete; implementation correctness and real-source compatibility
are not established for all workflows. Data availability and scientific validation
remain separate requirements:

- Scientific review of sample/feature/endpoint audits and additional full-source
  offline byte comparisons for the newly completed large source variants.
- Authorized original DKFZ/CPGEA sources if full raw-deposit coverage or an alternative
  to the documented public processed subsets is required.
- Verified patient-linked event/censoring times for cohorts whose public source only
  provides status or clinical covariates. Expression/PCA exports do not wait for these
  dates; survival fits remain unavailable until sufficient eligible data exist.
- Review of patient overlap, case-control/cohort sampling, tissue/treatment differences,
  expression normalization and targeted gene coverage before interpreting or pooling
  associations.

See [source and clinical blockers](remaining_datasets.md),
[clinical supplements](clinical_linkage.md), and the individual guides. A test pass or
structural checker pass is not evidence of clinical validity or a validated predictor.

## Review fixes and analysis policy

CancerMap now classifies the deposited `material_type` separately from organ `tissue`,
and excludes source-named xenografts. Its complete deposited-source build passed;
region sensitivity and clinical interpretation remain research decisions.
Clinical supplements require matching known endpoint definitions, origins and censoring rules;
unresolved deposited time definitions cannot be combined. Definitions for unmatched rows are preserved.
The twelve experimental cohorts export restrictive `analysis_policy` metadata.
Clinical record completeness is not analysis approval; downstream consumers must enforce
the recorded restrictions. This PR does not change the existing analysis commands.
See [review resolution](review_resolution.md) and [rebuild instructions](rebuilding.md).
