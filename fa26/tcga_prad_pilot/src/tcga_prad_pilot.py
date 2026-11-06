"""TCGA-PRAD pilot pipeline.

Downloads TCGA-PRAD gene expression and clinical data from cBioPortal, cleans both,
matches each primary-tumor RNA-seq sample to its patient's outcomes, runs QC, and saves
an analysis-ready dataset (three linked tables + manifest, dictionary, provenance, QC log).

Decisions applied (CampLab client meetings, Oct-Nov 2026):
  * Endpoint: progression-free interval (PFI) = cBioPortal PFS_STATUS / PFS_MONTHS
  * Matching key: first 12 characters of the TCGA barcode (patient); keep sample type 01
    (primary tumor), one sample per patient
  * Covariates: age, pathologic T stage, pathologic N stage, Gleason score
  * Genes with missing values in any sample are removed (all samples kept)
  * Samples more than z_threshold standard deviations from the mean on median expression or
    genes detected are excluded as QC outliers

Usage:
    python src/tcga_prad_pilot.py --out pilot_out [--module-file modules.csv] [--z-threshold 3]
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import io
import json
import os

import numpy as np
import pandas as pd
import requests

STUDY = "prad_tcga_pan_can_atlas_2018"   # TCGA PanCancer Atlas (cBioPortal)
LEGACY_STUDY = "prad_tcga"                # TCGA Firehose Legacy (used only for Gleason)
API = "https://www.cbioportal.org/api"
DATAHUB = "cBioPortal/datahub"
EXPR_FILE = f"public/{STUDY}/data_mrna_seq_v2_rsem.txt"
DETECTED_LOG2 = 1.0                       # a gene counts as "detected" above log2(RSEM + 1) = 1

DATA_DICTIONARY = [
    ("sample_id", "TCGA sample barcode (-01 = primary tumor)", "expression column names"),
    ("patient_id", "First 12 characters of the barcode; key to clinical data", "derived"),
    ("PFI_event", "Progression-free interval event: 1 = progression, 0 = censored", "PFS_STATUS"),
    ("PFI_months", "Months to progression or last follow-up", "PFS_MONTHS"),
    ("OS_event / OS_months", "Overall survival (secondary)", "OS_STATUS / OS_MONTHS"),
    ("DSS_event / DSS_months", "Disease-specific survival (secondary)", "DSS_STATUS / DSS_MONTHS"),
    ("DFI_event / DFI_months", "Disease-free interval (secondary; incomplete in PRAD)",
     "DFS_STATUS / DFS_MONTHS"),
    ("age", "Age at diagnosis (years)", "AGE"),
    ("path_T_stage / path_T_group", "Pathologic T stage (raw / grouped T2-T4)", "PATH_T_STAGE"),
    ("path_N_stage", "Pathologic N stage", "PATH_N_STAGE"),
    ("gleason_*", "Gleason fields from the Firehose Legacy study (blank tertiary = no tertiary pattern)",
     LEGACY_STUDY),
    ("expression values", "log2(RSEM normalized count + 1), one row per gene symbol", EXPR_FILE),
]


class QCLog:
    """Collects human-readable QC lines; printed as they are added and saved as qc_summary.md."""

    def __init__(self, verbose: bool = True):
        self.lines: list[str] = []
        self.verbose = verbose

    def note(self, msg: str) -> None:
        if self.verbose:
            print(msg)
        self.lines.append(msg)

    def text(self) -> str:
        return "# TCGA-PRAD pilot QC summary\n\n" + "\n".join(self.lines) + "\n"


# ------------------------------------------------------------------------- download
def api_get(path: str, **params):
    r = requests.get(f"{API}{path}", params=params, timeout=180)
    r.raise_for_status()
    return r.json()


def clinical_from_api(study: str, level: str) -> pd.DataFrame:
    """Clinical data as a wide table: one row per patient (PATIENT) or per sample (SAMPLE)."""
    raw = pd.DataFrame(api_get(f"/studies/{study}/clinical-data", clinicalDataType=level,
                               projection="SUMMARY", pageSize=100000))
    if level == "PATIENT":
        return raw.pivot(index="patientId", columns="clinicalAttributeId", values="value")
    wide = raw.pivot(index="sampleId", columns="clinicalAttributeId", values="value")
    wide["patientId"] = raw.drop_duplicates("sampleId").set_index("sampleId")["patientId"]
    return wide


def download_datahub_lfs(path: str) -> tuple[bytes, str]:
    """cBioPortal stores study files on GitHub with Git LFS. Try the media URL, then the LFS API."""
    media = f"https://media.githubusercontent.com/media/{DATAHUB}/master/{path}"
    try:
        r = requests.get(media, timeout=600)
        r.raise_for_status()
        if not r.content.startswith(b"version https://git-lfs"):
            return r.content, media
    except Exception as e:  # fall through to the LFS API
        print(f"  media URL failed ({e}); trying the Git LFS API")
    pointer = requests.get(f"https://raw.githubusercontent.com/{DATAHUB}/master/{path}", timeout=60).text
    oid = pointer.split("sha256:")[1].split()[0]
    size = int(pointer.split("size ")[1].split()[0])
    batch = requests.post(f"https://github.com/{DATAHUB}.git/info/lfs/objects/batch",
                          headers={"Accept": "application/vnd.git-lfs+json",
                                   "Content-Type": "application/vnd.git-lfs+json"},
                          json={"operation": "download", "transfers": ["basic"],
                                "objects": [{"oid": oid, "size": size}]}, timeout=60).json()
    href = batch["objects"][0]["actions"]["download"]["href"]
    r = requests.get(href, timeout=600)
    r.raise_for_status()
    return r.content, f"git-lfs sha256:{oid}"


def load_inputs() -> dict:
    """Download expression + clinical data and record where each came from."""
    prov = {"run_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "sources": {}}
    print("Downloading expression matrix (~77 MB) ...")
    expr_bytes, expr_src = download_datahub_lfs(EXPR_FILE)
    prov["sources"]["expression"] = {"study": STUDY, "file": EXPR_FILE, "url_or_oid": expr_src,
                                     "sha256": hashlib.sha256(expr_bytes).hexdigest(),
                                     "units": "RSEM normalized counts (batch-corrected PanCanAtlas)"}
    expr_raw = pd.read_csv(io.BytesIO(expr_bytes), sep="\t", low_memory=False)

    print("Downloading clinical data (cBioPortal API) ...")
    patients = clinical_from_api(STUDY, "PATIENT")
    samples = clinical_from_api(STUDY, "SAMPLE")
    prov["sources"]["clinical"] = {"study": STUDY, "api": f"{API}/studies/{STUDY}/clinical-data"}

    gleason = None
    try:
        legacy_p = clinical_from_api(LEGACY_STUDY, "PATIENT")
        legacy_s = clinical_from_api(LEGACY_STUDY, "SAMPLE")
        p_cols = [c for c in legacy_p.columns if "GLEASON" in c.upper()]
        s_cols = [c for c in legacy_s.columns if "GLEASON" in c.upper()]
        gleason = pd.concat([legacy_p[p_cols], legacy_s.groupby("patientId")[s_cols].first()], axis=1)
        prov["sources"]["gleason"] = {"study": LEGACY_STUDY, "columns": list(gleason.columns)}
    except Exception as e:
        print(f"  Gleason not available from {LEGACY_STUDY}: {e}")
    return {"expr_raw": expr_raw, "patients": patients, "samples": samples,
            "gleason": gleason, "provenance": prov}


# ------------------------------------------------------------------------- clean expression
def clean_expression(expr_raw: pd.DataFrame, log: QCLog):
    """Return (expr genes x samples in log2(RSEM + 1), feature table, list of removed genes)."""
    log.note("## Expression")
    df = expr_raw.copy()
    sample_cols = [c for c in df.columns if str(c).startswith("TCGA-")]
    log.note(f"- Raw matrix: {len(df):,} genes x {len(sample_cols)} samples")
    df["Hugo_Symbol"] = df["Hugo_Symbol"].astype(str).str.strip()
    bad_symbol = df["Hugo_Symbol"].isin(["", "nan", "?", "NA"])
    log.note(f"- Dropped {int(bad_symbol.sum())} rows with no gene symbol")
    df = df[~bad_symbol]

    values = df[sample_cols].apply(pd.to_numeric, errors="coerce").clip(lower=0)  # RSEM is never negative
    df = pd.concat([df[["Hugo_Symbol", "Entrez_Gene_Id"]], values], axis=1)
    df["_mean"] = values.mean(axis=1)
    n_dup = int(df["Hugo_Symbol"].duplicated(keep=False).sum())
    df = df.sort_values("_mean", ascending=False).drop_duplicates("Hugo_Symbol")
    log.note(f"- Duplicate gene symbols: {n_dup} rows involved; kept the most highly expressed row "
             f"per symbol ({len(df):,} genes)")

    features = df[["Hugo_Symbol", "Entrez_Gene_Id"]].set_index("Hugo_Symbol")
    expr = np.log2(df.set_index("Hugo_Symbol")[sample_cols] + 1)
    log.note("- Normalization: log2(RSEM normalized count + 1)")

    na_genes = list(expr.index[expr.isna().any(axis=1)])
    n_na_samples = int(expr.isna().any().sum())
    expr = expr.drop(index=na_genes)
    features = features.loc[expr.index]
    log.note(f"- Removed {len(na_genes):,} genes with missing values in {n_na_samples} samples "
             f"(PanCanAtlas batch-correction gaps); all samples kept")
    log.note(f"- Final: {len(expr):,} genes, no missing values")
    return expr, features, na_genes


# ------------------------------------------------------------------------- select samples
def select_samples(expr: pd.DataFrame, log: QCLog, z_threshold: float = 3.0) -> pd.DataFrame:
    """Sample manifest: barcode rule (type 01, one per patient) + QC outlier flags."""
    log.note("## Samples")
    man = pd.DataFrame({"sample_id": list(expr.columns)})
    man["patient_id"] = man["sample_id"].str[:12]
    man["sample_type_code"] = man["sample_id"].str[13:15]
    man["included"] = True
    man["reason"] = ""

    not_primary = man["sample_type_code"] != "01"
    man.loc[not_primary, ["included", "reason"]] = [False, "not primary tumor (type != 01)"]
    candidates = man[man["included"]].sort_values("sample_id")
    extra_idx = candidates.index[candidates.duplicated("patient_id")]
    man.loc[extra_idx, ["included", "reason"]] = [False, "extra sample for same patient"]
    dup_cols = man["sample_id"].duplicated()
    man.loc[dup_cols, ["included", "reason"]] = [False, "duplicate sample column"]
    log.note(f"- Expression samples: {len(man)}; non-primary removed: {int(not_primary.sum())}; "
             f"extra samples per patient removed: {len(extra_idx)}")

    kept = man.loc[man["included"], "sample_id"]
    med = expr[kept].median()
    detected = (expr[kept] > DETECTED_LOG2).sum()

    def z(s):
        return (s - s.mean()) / s.std(ddof=0)

    outlier = (z(med).abs() > z_threshold) | (z(detected).abs() > z_threshold)
    man["median_log2"] = man["sample_id"].map(med).round(3)
    man["genes_detected"] = man["sample_id"].map(detected)
    out_ids = outlier[outlier].index
    man.loc[man["sample_id"].isin(out_ids), ["included", "reason"]] = [
        False, f"QC outlier (|z|>{z_threshold:g})"]
    log.note(f"- QC outliers removed (median expression or genes detected, |z| > {z_threshold:g}): "
             f"{len(out_ids)}")
    return man


# ------------------------------------------------------------------------- clean clinical
def _event(series: pd.Series) -> pd.Series:
    """cBioPortal status strings look like '1:PROGRESSION' / '0:CENSORED'."""
    return series.astype(str).str[0].map({"1": 1, "0": 0})


def clean_clinical(patients: pd.DataFrame, samples: pd.DataFrame, gleason: pd.DataFrame | None,
                   log: QCLog) -> pd.DataFrame:
    """One row per patient: endpoints + covariates."""
    log.note("## Clinical")
    clin = pd.DataFrame(index=patients.index)
    clin["PFI_event"] = _event(patients["PFS_STATUS"])
    clin["PFI_months"] = pd.to_numeric(patients["PFS_MONTHS"], errors="coerce")
    for name, col in [("OS", "OS"), ("DSS", "DSS"), ("DFI", "DFS")]:
        if f"{col}_STATUS" in patients:
            clin[f"{name}_event"] = _event(patients[f"{col}_STATUS"])
            clin[f"{name}_months"] = pd.to_numeric(patients.get(f"{col}_MONTHS"), errors="coerce")
    clin["age"] = pd.to_numeric(patients.get("AGE"), errors="coerce")
    for src, dst in [("PATH_T_STAGE", "path_T_stage"), ("PATH_N_STAGE", "path_N_stage")]:
        col = patients.get(src)
        if col is None and src in samples:
            col = samples.groupby("patientId")[src].first()
        clin[dst] = col
    clin["path_T_group"] = clin["path_T_stage"].astype(str).str.extract(r"(T\d)", expand=False)
    if gleason is not None:
        for c in gleason.columns:
            clin[c.lower()] = gleason[c].reindex(clin.index)
    bad_time = clin["PFI_months"] < 0
    clin.loc[bad_time, "PFI_months"] = np.nan
    log.note(f"- Patients: {len(clin)}; negative PFI times set to missing: {int(bad_time.sum())}")
    return clin


# ------------------------------------------------------------------------- match
def match(manifest: pd.DataFrame, clin: pd.DataFrame, expr: pd.DataFrame, log: QCLog):
    """Join included samples to patient outcomes. Returns (clinical, expression, manifest)."""
    log.note("## Matching")
    man = manifest.copy()
    no_clin = man["included"] & ~man["patient_id"].isin(clin.index)
    man.loc[no_clin, ["included", "reason"]] = [False, "no clinical record"]

    final = man[man["included"]].set_index("sample_id")
    clinical = final[["patient_id"]].join(clin, on="patient_id")
    expression = expr[final.index]
    assert clinical["patient_id"].is_unique, "a patient appears more than once"
    assert list(expression.columns) == list(clinical.index), "expression and clinical are out of order"

    usable = clinical["PFI_event"].notna() & clinical["PFI_months"].notna()
    n_out = len(set(clin.index) - set(clinical["patient_id"]))
    n_qc = int((man["reason"].str.startswith("QC") & man["patient_id"].isin(clin.index)).sum())
    log.note(f"- Final matched samples: {len(clinical)} (1 per patient, all primary tumor)")
    log.note(f"- Patients not in final dataset: {n_out} ({n_out - n_qc} with no RNA-seq, "
             f"{n_qc} QC exclusions)")
    log.note(f"- Usable for PFI: {int(usable.sum())}; events: {int(clinical.loc[usable, 'PFI_event'].sum())}; "
             f"median follow-up {clinical.loc[usable, 'PFI_months'].median():.1f} months")
    miss = (clinical.isna().mean() * 100).round(1)
    log.note("- Missing (%): " + ", ".join(f"{k} {v}" for k, v in miss.items() if v > 0))
    return clinical, expression, man


def annotate_features(features: pd.DataFrame, expression: pd.DataFrame) -> pd.DataFrame:
    f = features.loc[expression.index].copy()
    f["mean_log2"] = expression.mean(axis=1).round(3)
    f["pct_samples_expressed"] = ((expression > DETECTED_LOG2).mean(axis=1) * 100).round(1)
    return f


# ------------------------------------------------------------------------- module coverage
def module_coverage(module_file: str, genes, removed_genes, log: QCLog) -> pd.DataFrame:
    """Coverage of each module's genes. File: CSV/TSV with a gene column and a module column."""
    mods = pd.read_csv(module_file, sep=None, engine="python")
    gcol = next(c for c in mods.columns
                if c.lower() in ("gene", "genes", "symbol", "gene_symbol", "hugo_symbol"))
    mcol = next(c for c in mods.columns if c != gcol)
    sym = mods[gcol].astype(str).str.strip()
    mods["in_tcga"] = sym.isin(set(genes))
    mods["removed_missing_values"] = sym.isin(set(removed_genes))
    cov = mods.groupby(mcol).agg(genes=("in_tcga", "size"), found=("in_tcga", "sum"),
                                 lost_to_missing_values=("removed_missing_values", "sum"))
    cov["pct_found"] = (100 * cov["found"] / cov["genes"]).round(1)
    log.note("## Module coverage\n" + cov.to_string())
    return cov


# ------------------------------------------------------------------------- save
def save_outputs(out_dir, expression, clinical, features, manifest, provenance, removed_genes,
                 log: QCLog, coverage: pd.DataFrame | None = None) -> None:
    os.makedirs(out_dir, exist_ok=True)
    expression.round(4).to_csv(f"{out_dir}/expression_log2_genes_x_samples.csv.gz")
    clinical.to_csv(f"{out_dir}/clinical_by_sample.csv")
    features.to_csv(f"{out_dir}/feature_annotation.csv")
    manifest.to_csv(f"{out_dir}/sample_manifest.csv", index=False)
    pd.Series(removed_genes, name="gene").to_csv(f"{out_dir}/genes_removed_missing_values.csv", index=False)
    pd.DataFrame(DATA_DICTIONARY, columns=["field", "description", "source"]).to_csv(
        f"{out_dir}/data_dictionary.csv", index=False)
    if coverage is not None:
        coverage.to_csv(f"{out_dir}/module_coverage.csv")
    with open(f"{out_dir}/provenance.json", "w") as fh:
        json.dump(provenance, fh, indent=2)
    log.note("## Output")
    try:
        import anndata as ad
        obs = clinical.copy()
        text_cols = obs.select_dtypes(exclude="number").columns
        obs[text_cols] = obs[text_cols].astype(str)
        var = features.copy()
        var["Entrez_Gene_Id"] = var["Entrez_Gene_Id"].astype(str)
        adata = ad.AnnData(X=expression.T.values.astype("float32"), obs=obs, var=var)
        adata.uns["provenance"] = json.dumps(provenance)
        adata.write_h5ad(f"{out_dir}/tcga_prad_pilot.h5ad")
        log.note(f"- AnnData: {adata.n_obs} samples x {adata.n_vars:,} genes")
    except ImportError:
        log.note("- AnnData not saved (install `anndata`); CSV tables are complete")
    with open(f"{out_dir}/qc_summary.md", "w") as fh:
        fh.write(log.text())


# ------------------------------------------------------------------------- run everything
def run(out_dir: str = "pilot_out", module_file: str | None = None, z_threshold: float = 3.0,
        verbose: bool = True) -> dict:
    log = QCLog(verbose)
    inputs = load_inputs()
    expr, features, removed = clean_expression(inputs["expr_raw"], log)
    manifest = select_samples(expr, log, z_threshold)
    clin = clean_clinical(inputs["patients"], inputs["samples"], inputs["gleason"], log)
    clinical, expression, manifest = match(manifest, clin, expr, log)
    features = annotate_features(features, expression)
    coverage = module_coverage(module_file, expression.index, removed, log) if module_file else None
    save_outputs(out_dir, expression, clinical, features, manifest, inputs["provenance"], removed, log,
                 coverage)
    print(f"\nDone. Outputs in ./{out_dir}/")
    return {"expression": expression, "clinical": clinical, "features": features,
            "manifest": manifest, "removed_genes": removed, "coverage": coverage, "log": log}


def main() -> None:
    ap = argparse.ArgumentParser(description="TCGA-PRAD pilot pipeline")
    ap.add_argument("--out", default="pilot_out", help="output folder (default: pilot_out)")
    ap.add_argument("--module-file", default=None, help="CSV/TSV with gene and module columns")
    ap.add_argument("--z-threshold", type=float, default=3.0, help="QC outlier cutoff (default: 3)")
    args, _ = ap.parse_known_args()  # tolerate extra args passed by Jupyter/Colab
    run(args.out, args.module_file, args.z_threshold)


if __name__ == "__main__":
    main()
