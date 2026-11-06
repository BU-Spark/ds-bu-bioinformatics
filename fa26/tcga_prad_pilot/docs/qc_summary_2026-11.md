# TCGA-PRAD pilot QC summary (November 2026 run)

## Expression
- Raw matrix: 20,531 genes x 493 samples
- Dropped 0 rows with no gene symbol
- Duplicate gene symbols: 27 rows involved; kept the most highly expressed row per symbol (20,512 genes)
- Normalization: log2(RSEM normalized count + 1)
- Removed 1,439 genes with missing values in 33 samples (PanCanAtlas batch-correction gaps); all samples kept
- Final: 19,073 genes, no missing values

## Samples
- Expression samples: 493; non-primary removed: 0; extra samples per patient removed: 0
- QC outliers removed (median expression or genes detected, |z| > 3): 9 (includes 2 PFI events; event rate similar to the cohort)

## Clinical
- Patients: 494; negative PFI times set to missing: 0

## Matching
- Final matched samples: 484 (1 per patient, all primary tumor)
- Patients not in final dataset: 10 (1 with no RNA-seq, 9 QC exclusions)
- Usable for PFI: 484; events: 91; median follow-up 25.8 months
- Missing (%): DSS_event 0.4, DFI_event 31.8, DFI_months 31.8, path_T_stage 1.4, path_N_stage 14.7, path_T_group 1.4, gleason_pattern_tertiary 69.2 (blank = no tertiary pattern)

## Output
- AnnData: 484 samples x 19,073 genes
