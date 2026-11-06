"""Feature-preserving GSE16560 (Sboner / Swedish Watchful-Waiting) pipeline.

GSE16560-specific counterpart of Mansur Orazkhan's Taylor/GSE21034 workflow
(BU-Spark/ds-bu-bioinformatics, branch mansur_orazkhan). Same design principles:

* checksummed raw inputs recorded in a manifest (reused/verified on re-runs)
* strict parsing of the deposited Series Matrix; no extra normalization,
  imputation, filtering or probe averaging
* every GPL5474 probe stays a separate feature; gene IDs/symbols may repeat
* audit tables for samples and features, explicit nullable clinical fields
* CSV + AnnData (.h5ad) exports validated by reloading, plus metadata.json,
  validation.json and checksums.sha256
* existing differing outputs are protected, never overwritten

Dataset: Illumina DASL Human 6k (GPL5474), FFPE tumor, 281 untreated localized
PCa cases (TURP at diagnosis), extreme-groups design (165 lethal / 116 indolent).

Differences from Taylor: the platform provides GeneID/Symbol/RefSeq directly, so
no MyGene lookup is needed (annotation is validated, not re-derived); clinical
data come from the GEO sample characteristics (no cBioPortal).

Usage: uv run --locked sboner-gse16560 [--data-dir data] [--offline]
"""
from __future__ import annotations
import argparse
import csv
import gzip
import io
import json
import re
import tempfile
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from .anndata_export import write_anndata
from .expression import write_csv
from .sources import SourceCache, sha256, write_json
from .taylor_pipeline import check_destination, publish

ACC, PLATFORM = "GSE16560", "GPL5474"
MATRIX_NAME, GPL_NAME = f"{ACC}_series_matrix.txt.gz", f"{PLATFORM}_family.soft.gz"
MATRIX_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/{ACC}/matrix/{MATRIX_NAME}"
GPL_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL5nnn/{PLATFORM}/soft/{GPL_NAME}"
N_TRAIN = 186  # series overall_design: first 186 samples = training, remaining 95 = validation
EXPECTED_SAMPLES = 281
OUTCOME_CODES = {"lethal": 1, "indolent": 0}  # extreme-groups label -> PCSM indicator
VITAL_CODES = {"dead": 1, "alive": 0}
REQUIRED_CHARS = ["gleason", "major.gleason", "minor.gleason", "batch", "extreme", "cancer.percent",
                  "age", "diag.yr", "status.all", "fup.month", "fusion"]
SNAPSHOT_DEFAULT = "sboner-v1"


# --------------------------------------------------------------------------- expression
def read_geo(path: Path):
    """Strictly parse the Series Matrix -> (sample metadata, probes x GSM expression, processing text)."""
    rows, chars, series, begin = {}, [], None, False
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0] == "!series_matrix_table_begin":
                begin = True
                break
            if row[0] == "!Series_geo_accession":
                series = row[1:]
            elif row[0] == "!Sample_characteristics_ch1":
                chars.append(row[1:])
            elif row[0] in {"!Sample_geo_accession", "!Sample_title", "!Sample_platform_id", "!Sample_data_processing"}:
                if row[0] in rows:
                    raise ValueError(f"Duplicate GEO header: {row[0]}")
                rows[row[0]] = row[1:]
        if series != [ACC]:
            raise ValueError(f"Expected GEO series {ACC}")
        if not begin or len(rows) != 4:
            raise ValueError("Missing GEO header or table begin")
        acc = rows["!Sample_geo_accession"]
        n = len(acc)
        if not n or len(set(acc)) != n or any(not x for x in acc):
            raise ValueError("Duplicate/empty GEO accessions")
        if any(len(v) != n for v in [*rows.values(), *chars]):
            raise ValueError("Inconsistent GEO metadata lengths")
        if set(rows["!Sample_platform_id"]) != {PLATFORM}:
            raise ValueError(f"Expected expression platform {PLATFORM}")
        records = []
        for i in range(n):
            rec = {"geo_accession": acc[i], "title": rows["!Sample_title"][i], "platform_id": PLATFORM}
            for values in chars:
                if not values[i].strip():
                    continue
                if ":" not in values[i]:
                    raise ValueError("Malformed GEO characteristic")
                key, value = values[i].split(":", 1)
                key = key.strip()
                if key in rec:
                    raise ValueError(f"Duplicate characteristic: {key}")
                rec[key] = value.strip()
            records.append(rec)
        meta = pd.DataFrame(records)
        if meta.title.duplicated().any():
            raise ValueError("Sample titles (patient IDs) are not unique")
        header = next(reader, [])
        if header != ["ID_REF", *acc]:
            raise ValueError("Matrix columns differ from GEO accessions/order")
        features, values, end = [], [], False
        for row in reader:
            if row == ["!series_matrix_table_end"]:
                end = True
                break
            if len(row) != n + 1 or not row[0]:
                raise ValueError("Malformed expression row")
            features.append(row[0])
            values.append(row[1:])
        if not end or not features or len(set(features)) != len(features):
            raise ValueError("Truncated/empty matrix or duplicate feature IDs")
    expr = pd.DataFrame(np.asarray(values, dtype=float), index=pd.Index(features, name="feature_id"), columns=acc)
    if not np.isfinite(expr.to_numpy()).all():
        raise ValueError("Missing/non-finite expression; no implicit imputation")
    return meta, expr, sorted(set(rows["!Sample_data_processing"]))


# --------------------------------------------------------------------------- annotation
def parse_platform(path: Path) -> pd.DataFrame:
    """GPL5474 table -> one row per probe with platform-supplied GeneID/Symbol/RefSeq."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        text = f.read()
    if f"^PLATFORM = {PLATFORM}" not in text or "!platform_table_end" not in text:
        raise ValueError(f"Wrong/truncated {PLATFORM} SOFT record")
    table = text.split("!platform_table_begin\n", 1)[1].split("!platform_table_end", 1)[0]
    ann = pd.read_csv(io.StringIO(table), sep="\t", dtype=str, keep_default_na=False)
    if not {"ID", "GeneID", "Symbol", "REFSEQ", "GB_ACC"}.issubset(ann.columns):
        raise ValueError("Unexpected GPL5474 annotation columns")
    ann = ann.rename(columns={"ID": "feature_id", "GeneID": "gene_id", "Symbol": "gene_symbol",
                              "REFSEQ": "refseq_all", "GB_ACC": "refseq_original"})
    if ann.feature_id.duplicated().any() or ann.feature_id.eq("").any():
        raise ValueError("Duplicate/empty GPL5474 probe IDs")
    ann["gene_symbol"] = ann.gene_symbol.str.strip()
    ann["refseq_query"] = ann.refseq_original.str.replace(r"\.\d+$", "", regex=True)
    return ann[["feature_id", "gene_id", "gene_symbol", "refseq_original", "refseq_query", "refseq_all"]]


def resolve_mapping(platform: pd.DataFrame, expression_ids) -> pd.DataFrame:
    """Audit every probe; accept only a valid, internally consistent platform annotation.

    Accepted: numeric Entrez ID, non-empty symbol, RefSeq mRNA accession, and a
    one-to-one Entrez<->symbol identity across the platform. Everything else is
    excluded with a reason (never silently dropped, never selected by score).
    """
    if set(expression_ids) != set(platform.feature_id):
        raise ValueError("Expression/GPL feature sets differ; no silent inner join")
    m = platform.copy()
    valid_id = m.gene_id.str.fullmatch(r"\d+")
    valid_sym = m.gene_symbol.ne("")
    valid_ref = m.refseq_query.str.fullmatch(r"(?:NM|NR|XM|XR)_\d+")
    in_list = pd.Series([q in {re.sub(r"\.\d+$", "", x.strip()) for x in a.split(";")}
               for q, a in zip(m.refseq_query, m.refseq_all)], index=m.index)
    status = pd.Series("mapped", index=m.index)
    status[~in_list] = "refseq_not_in_platform_list"
    status[~valid_ref] = "invalid_refseq"
    status[~valid_sym] = "no_symbol"
    status[~valid_id] = "no_entrez_id"
    ok = status.eq("mapped")
    pairs = m.loc[ok, ["gene_id", "gene_symbol"]].drop_duplicates()
    bad_sym = set(pairs.loc[pairs.gene_symbol.duplicated(False), "gene_symbol"])
    bad_id = set(pairs.loc[pairs.gene_id.duplicated(False), "gene_id"])
    status[ok & (m.gene_symbol.isin(bad_sym) | m.gene_id.isin(bad_id))] = "gene_symbol_identity_conflict"
    m["mapping_status"], m["selected"] = status, status.eq("mapped")
    m["n_probes_same_gene"] = m.groupby("gene_id")["feature_id"].transform("count").where(m.selected, 0)
    m["mapping_source"] = f"{PLATFORM} platform annotation (GeneID/Symbol/GB_ACC); not re-derived"
    return m


def select_features(expression: pd.DataFrame, mapping: pd.DataFrame):
    chosen = mapping[mapping.selected].set_index("feature_id").loc[lambda d: d.index.isin(expression.index)]
    if chosen.empty:
        raise ValueError("No mapped features")
    chosen = chosen.loc[expression.index[expression.index.isin(chosen.index)]]  # keep source order
    selected = expression.loc[chosen.index].copy()
    genes = chosen[["gene_id", "gene_symbol", "refseq_original", "refseq_query", "n_probes_same_gene"]].reset_index()
    return selected, genes


# --------------------------------------------------------------------------- clinical
def _num(series: pd.Series, name: str) -> pd.Series:
    out = pd.to_numeric(series.replace({"NA": np.nan, "": np.nan}), errors="coerce")
    bad = series.notna() & ~series.isin(["NA", ""]) & out.isna()
    if bad.any():
        raise ValueError(f"Non-numeric values in {name}: {sorted(series[bad].unique())[:5]}")
    return out


def build_clinical(meta: pd.DataFrame):
    missing = [c for c in REQUIRED_CHARS if c not in meta]
    if missing:
        raise ValueError(f"Missing GEO characteristics: {missing}")
    outcome = meta["extreme"].str.lower()
    vital = meta["status.all"].str.lower()
    if not set(outcome) <= set(OUTCOME_CODES) or not set(vital) <= set(VITAL_CODES):
        raise ValueError(f"Unknown outcome/vital codes: {set(outcome) - set(OUTCOME_CODES)} {set(vital) - set(VITAL_CODES)}")
    fup = _num(meta["fup.month"], "fup.month")
    if (fup < 0).any() or not np.isfinite(fup.dropna()).all():
        raise ValueError("Negative/non-finite follow-up time")
    c = pd.DataFrame({
        "sample_id": meta.geo_accession, "geo_accession": meta.geo_accession, "patient_id": meta.title,
        "platform_id": meta.platform_id,
        "age": _num(meta["age"], "age"), "gleason_score": _num(meta["gleason"], "gleason"),
        "gleason_primary": _num(meta["major.gleason"], "major.gleason"),
        "gleason_secondary": _num(meta["minor.gleason"], "minor.gleason"),
        "diagnosis_year": _num(meta["diag.yr"], "diag.yr"),
        "cancer_percent": _num(meta["cancer.percent"], "cancer.percent"),
        "batch": meta["batch"], "tmprss2_erg_fusion": _num(meta["fusion"], "fusion"),
        "outcome_group": outcome, "vital_status": vital, "followup_months": fup,
    })
    c["pcsm_event"] = outcome.map(OUTCOME_CODES).astype(int)
    c["os_event"] = vital.map(VITAL_CODES).astype(int)
    order = np.arange(len(c))
    c["split"] = np.where(order < N_TRAIN, "train", "validation")  # relies on deposited sample order
    c["survival_eligible"] = c.followup_months.notna() & (c.followup_months > 0)
    c["survival_ineligible_reason"] = np.where(c.followup_months.isna(), "missing_followup",
                                       np.where(c.followup_months <= 0, "zero_followup", ""))
    c["stage"] = "T1-T2 N0 (cohort-level; per-patient stage not in GEO)"
    c["treatment"] = "none (TURP at diagnosis; watchful waiting)"
    c["tissue"] = "FFPE tumor, primary"
    return c.reset_index(drop=True)


def build_audit(meta: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({"geo_accession": meta.geo_accession, "patient_id": meta.title,
                         "platform_id": meta.platform_id, "selected": True,
                         "reason": "cohort is 100% primary FFPE tumor, one sample per patient",
                         "outcome_group": clinical.outcome_group})


# --------------------------------------------------------------------------- validation
def validate(expression, clinical, samples, genes, mapping, audit) -> dict:
    if list(expression.columns) != list(clinical.sample_id) or list(clinical.sample_id) != list(samples.sample_id):
        raise ValueError("Expression/clinical/sample order mismatch")
    if list(expression.index) != list(genes.feature_id):
        raise ValueError("Expression/gene annotation order mismatch")
    if not expression.index.is_unique or not clinical.sample_id.is_unique or not clinical.patient_id.is_unique:
        raise ValueError("Non-unique feature/sample/patient IDs")
    if not np.isfinite(expression.to_numpy()).all():
        raise ValueError("Non-finite expression")
    if len(clinical) != EXPECTED_SAMPLES:
        raise ValueError(f"Expected {EXPECTED_SAMPLES} samples, got {len(clinical)}")
    counts = {"samples": len(clinical), "patients": int(clinical.patient_id.nunique()),
              "source_features": len(mapping), "retained_feature_rows": len(expression),
              "unique_gene_ids": int(genes.gene_id.nunique()), "unique_gene_symbols": int(genes.gene_symbol.nunique()),
              "features_sharing_gene_id": int((genes.n_probes_same_gene > 1).sum()),
              "pcsm_events": int(clinical.pcsm_event.sum()), "indolent": int((1 - clinical.pcsm_event).sum()),
              "deaths_any_cause": int(clinical.os_event.sum()),
              "train": int((clinical.split == "train").sum()), "validation": int((clinical.split == "validation").sum())}
    return {"counts": counts,
            "excluded_features_by_status": mapping.loc[~mapping.selected, "mapping_status"].value_counts().to_dict(),
            "expression_range": [float(expression.to_numpy().min()), float(expression.to_numpy().max())],
            "missing_clinical_values": {k: int(v) for k, v in clinical.isna().sum().items() if v},
            "sample_order_aligned": True, "feature_order_aligned": True,
            "outcome_group_vs_pcsm_event_consistent":
                bool((clinical.outcome_group.eq("lethal") == clinical.pcsm_event.eq(1)).all())}


def validate_export(stage: Path, expression: pd.DataFrame) -> dict:
    back = pd.read_csv(stage / "expression.csv.gz", index_col=0, float_precision="round_trip")
    if list(back.index) != list(expression.index) or list(back.columns) != list(expression.columns):
        raise ValueError("CSV round trip changed order")
    if not np.array_equal(back.to_numpy(), expression.to_numpy()):
        raise ValueError("CSV round trip changed expression values")
    clin = pd.read_csv(stage / "clinical.csv")
    genes = pd.read_csv(stage / "genes.csv", dtype=str)
    return {"counts": {"samples": len(clin), "retained_feature_rows": len(genes)}}


def build(cache: SourceCache, output_root: Path, interim_root: Path, anndata_root: Path,
          snapshot: str = SNAPSHOT_DEFAULT) -> dict:
    cache.verify_all()
    matrix_path = cache.fetch(MATRIX_NAME, MATRIX_URL)
    gpl_path = cache.fetch(GPL_NAME, GPL_URL)
    meta, features, processing = read_geo(matrix_path)
    platform = parse_platform(gpl_path)
    clinical = build_clinical(meta)
    audit = build_audit(meta, clinical)
    samples = clinical[["sample_id", "geo_accession", "patient_id", "platform_id"]].copy()
    features.columns = pd.Index(clinical.sample_id, name=None)  # GSM IDs are the sample IDs
    mapping = resolve_mapping(platform, features.index)
    expression, genes = select_features(features, mapping)
    checks = validate(expression, clinical, samples, genes, mapping, audit)
    checks["expression_scale_verified"] = False  # log2 stated by GEO text, not independently verified
    checks["pcsm_time_origin_verified"] = False

    metadata_out = {
        "schema_version": 1, "study": "Sboner / Swedish Watchful Waiting", "geo_series": ACC, "platform": PLATFORM,
        "snapshot": snapshot, "source_manifest_sha256": sha256(cache.manifest_path),
        "expression": {"orientation": "features_by_samples", "feature_key": "original GPL5474 probe ID; gene IDs/symbols may repeat",
                       "source_processing": processing,
                       "scale": "deposited signal; GEO text states log2, median-centered, qspline-normalized within DAP",
                       "additional_transformations": "none", "aggregation": "none; each mapped probe retained separately"},
        "mapping": {"source": f"{PLATFORM} platform table (GeneID, Symbol, GB_ACC)",
                    "policy": "Numeric Entrez ID, non-empty symbol, valid RefSeq accession listed in platform REFSEQ field, one-to-one Entrez/symbol identity; exclusions recorded with reasons",
                    "symbol_validation": "platform labels; not claimed HGNC-approved"},
        "outcome": {"endpoint": "PCSM (prostate-cancer-specific mortality)",
                    "event_definition": "pcsm_event = 1 for GEO 'extreme: Lethal', 0 for 'Indolent'",
                    "design_warning": "extreme-groups case-control sampling: lethal/indolent proportions are set by design, not cohort incidence; event rates are not population estimates and standard Cox baseline hazards are biased",
                    "time": "followup_months from GEO fup.month; time origin (diagnosis vs TURP) not stated in GEO",
                    "overall_survival": "os_event from status.all (all-cause death)"},
        "sample_policy": "All 281 deposited samples (single FFPE primary tumor per patient); no exclusions",
        "split": f"first {N_TRAIN} deposited samples = train, remaining = validation (per series overall_design)",
        "overlap_warning": "Overlaps the Swedish WW population in Setlur (GSE8402); not independent",
        "limitations": ["Targeted ~6k-gene DASL panel, not genome-wide", "Per-patient stage not deposited",
                        "Log scale taken from GEO text; no raw intensity QC", "PCSM time origin unverified",
                        "Extreme-groups design"],
        "code_sha256": {p.name: sha256(p) for p in sorted(Path(__file__).parent.glob("*.py"))},
    }
    metadata_out["output_schemas"] = {name: {c: {"dtype": str(f[c].dtype), "nullable": bool(f[c].isna().any())} for c in f.columns}
                                      for name, f in [("clinical.csv", clinical), ("samples.csv", samples), ("genes.csv", genes),
                                                      ("feature_mapping.csv.gz", mapping), ("sample_audit.csv", audit)]}
    for root in (output_root, interim_root, anndata_root):
        root.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".g16560-", dir=output_root))
    interim_stage = Path(tempfile.mkdtemp(prefix=".g16560-", dir=interim_root))
    anndata_stage = Path(tempfile.mkdtemp(prefix=".g16560-", dir=anndata_root))
    try:
        write_csv(expression, stage / "expression.csv.gz", index=True)
        for frame, name in [(clinical, "clinical.csv"), (samples, "samples.csv"), (genes, "genes.csv"),
                            (mapping, "feature_mapping.csv.gz"), (audit, "sample_audit.csv")]:
            write_csv(frame, stage / name)
        write_csv(features, interim_stage / "feature_expression.csv.gz", index=True)
        if validate_export(stage, expression)["counts"] != {"samples": checks["counts"]["samples"],
                                                              "retained_feature_rows": checks["counts"]["retained_feature_rows"]}:
            raise ValueError("Export round-trip changed counts")
        checks["export_round_trip_valid"] = True
        write_anndata(expression, clinical, genes, metadata_out, anndata_stage / "expression.h5ad")
        checks["anndata_round_trip_valid"] = True
        write_json(anndata_stage / "metadata.json", metadata_out)
        (anndata_stage / "checksums.sha256").write_text("".join(f"{sha256(p)}  {p.name}\n" for p in sorted(anndata_stage.iterdir())))
        write_json(stage / "metadata.json", metadata_out)
        write_json(stage / "validation.json", checks)
        (stage / "checksums.sha256").write_text("".join(f"{sha256(p)}  {p.name}\n" for p in sorted(stage.iterdir())))
        for staged, root in [(interim_stage, interim_root), (stage, output_root), (anndata_stage, anndata_root)]:
            check_destination(staged, root / snapshot)
        publish(interim_stage, interim_root / snapshot)
        publish(stage, output_root / snapshot)
        publish(anndata_stage, anndata_root / snapshot)
    finally:
        for s in (stage, interim_stage, anndata_stage):
            shutil.rmtree(s, ignore_errors=True)
    return checks


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", type=Path, default=Path("data"))
    p.add_argument("--raw-root", type=Path, help="raw snapshot dir (default <data-dir>/raw)")
    p.add_argument("--snapshot", default=SNAPSHOT_DEFAULT, help="output snapshot name")
    p.add_argument("--offline", action="store_true", help="forbid HTTP; require cached, checksummed inputs")
    p.add_argument("--output-root", type=Path)
    p.add_argument("--interim-root", type=Path)
    p.add_argument("--anndata-root", type=Path)
    a = p.parse_args(argv)
    if not a.snapshot or Path(a.snapshot).name != a.snapshot or a.snapshot in {".", ".."}:
        p.error("snapshot must be a single directory name")
    raw = a.raw_root or a.data_dir / "raw" / a.snapshot
    out = a.output_root or a.data_dir / "processed"
    interim = a.interim_root or a.data_dir / "interim"
    ann = a.anndata_root or a.data_dir / "anndata"
    paths = [r.resolve() for r in (raw, out, interim, ann)]
    if any(x.is_relative_to(y) or y.is_relative_to(x) for i, x in enumerate(paths) for y in paths[i + 1:]):
        p.error("raw, output, interim and AnnData directories must not overlap")
    cache = SourceCache(raw, offline=a.offline)
    if not cache.manifest['sources']:  # SourceCache defaults to the Taylor accessions
        cache.manifest['accessions'] = [ACC, PLATFORM]
    print(json.dumps(build(cache, out, interim, ann, a.snapshot), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
