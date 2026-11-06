# TCGA-PRAD pilot pipeline

Pilot data pipeline for the **CampLab Modular Cancer Cell Analysis** project (BU Spark!, DS 701, Fall 2026).
It downloads TCGA prostate adenocarcinoma (TCGA-PRAD) gene expression and clinical data, cleans both,
matches each primary-tumor RNA-seq sample to its patient's outcomes, runs quality control, and saves an
analysis-ready dataset. It is the template for harmonizing the other validation cohorts.

## Quick start

All commands run from this folder (`fa26/tcga_prad_pilot/`):

```bash
cd fa26/tcga_prad_pilot
pip install -r requirements.txt
python src/tcga_prad_pilot.py                       # outputs go to ./pilot_out/
python src/tcga_prad_pilot.py --module-file modules.csv   # also report module gene coverage
pytest -q                                           # offline tests on synthetic data
```

Or open `notebooks/tcga_prad_pilot_pipeline.ipynb`, which runs the same steps one at a time so each can be
inspected. In Google Colab: `!git clone -b patrick_mulcahy https://github.com/BU-Spark/ds-bu-bioinformatics.git`, `%cd ds-bu-bioinformatics/fa26/tcga_prad_pilot/notebooks`, then *Runtime → Run all*.

A full run takes about 3–5 minutes and needs internet access (cBioPortal and GitHub).

## Data sources

| Data | Source |
|---|---|
| Expression | cBioPortal study `prad_tcga_pan_can_atlas_2018`, file `data_mrna_seq_v2_rsem.txt` (RSEM normalized counts, PanCanAtlas batch-corrected), read from the cBioPortal datahub on GitHub |
| Clinical outcomes and covariates | cBioPortal REST API, same study |
| Gleason score | cBioPortal study `prad_tcga` (TCGA Firehose Legacy) |

Each run records the exact source and a SHA-256 fingerprint of the expression file in `provenance.json`.

## Pipeline steps and decisions

| Step | What it does |
|---|---|
| 1. Download | Expression matrix + patient- and sample-level clinical data |
| 2. Clean expression | Drop rows without a gene symbol; keep the most highly expressed row per duplicate symbol; transform to log2(RSEM + 1); **remove genes with missing values** (PanCanAtlas batch-correction gaps) so all samples are kept |
| 3. Select samples | **Barcode rule:** characters 1–12 = patient, 14–15 = sample type; keep `01` (primary tumor), one sample per patient. **QC:** exclude samples more than 3 SD from the mean on median expression or genes detected |
| 4. Clean clinical | **Endpoint: progression-free interval (PFI)** from `PFS_STATUS` / `PFS_MONTHS`. Secondary: OS, DSS, DFI. Covariates: age, pathologic T and N stage, Gleason |
| 5. Match | Join samples to patients on the 12-character ID; stops with an error on duplicate patients or misaligned tables |
| 6. Module coverage (optional) | Share of each module's genes measured, and how many were lost to the missing-value filter |
| 7. Save | Writes the outputs below |

Endpoint, matching rule and covariates follow guidance from the client (CampLab) meetings.

## Outputs (`pilot_out/`, not committed to git)

| File | Contents |
|---|---|
| `expression_log2_genes_x_samples.csv.gz` | Expression matrix, genes × samples |
| `clinical_by_sample.csv` | Endpoints and covariates, one row per sample |
| `feature_annotation.csv` | Gene symbol, Entrez ID, mean expression, % samples expressed |
| `tcga_prad_pilot.h5ad` | The three tables above combined in one AnnData object |
| `sample_manifest.csv` | Every sample, whether it was kept, and why not |
| `genes_removed_missing_values.csv` | Genes dropped by the missing-value filter |
| `data_dictionary.csv` | Definition and source of each field |
| `provenance.json` | Data sources, file fingerprint and run time |
| `qc_summary.md` | Step-by-step counts from the run |

## Results (November 2026 run)

- **484 samples × 19,073 genes**, one primary tumor per patient, no missing expression values
- **91 PFI events**, median follow-up 25.8 months
- 1,439 genes removed for missing values in 33 samples; 9 samples (1.8%) excluded as QC outliers, with an
  event rate similar to the rest of the cohort

Full counts: [`docs/qc_summary_2026-11.md`](docs/qc_summary_2026-11.md).

## Repository layout

```
src/tcga_prad_pilot.py               pipeline (importable functions + command-line entry point)
notebooks/tcga_prad_pilot_pipeline.ipynb   step-by-step walkthrough
tests/test_pipeline.py               offline tests on synthetic data
docs/qc_summary_2026-11.md           QC summary from the November 2026 run
```

## Open items

- Module coverage check once the module gene lists are available (`--module-file`)
- Client confirmation of the missing-gene filter and QC exclusions
- Extend the same structure to the other Tier 1 cohorts (Belfast, Sboner, Erho), which need probe-to-gene mapping
