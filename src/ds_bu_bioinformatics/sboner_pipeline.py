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
from .provenance import code_hashes
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
SNAPSHOT_DEFAULT = "sboner-v2"


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
MISSING_TOKENS = {"", "NA", "N/A", "NAN", "NULL", "[NOT AVAILABLE]"}
# endpoint prefix -> (source status column, source time column, status codes); same column
# contract as the other cohorts: <e>_event, <e>_months, <e>_eligible, <e>_eligibility_reason, <e>_time_unit
ENDPOINTS = {"pcsm": ("outcome_group", "followup_months", OUTCOME_CODES),
             "os": ("vital_status", "followup_months", VITAL_CODES)}
PRIMARY_ENDPOINT = "pcsm"


def _is_missing(value) -> bool:
    return pd.isna(value) or str(value).strip().upper() in MISSING_TOKENS


def parse_number(series: pd.Series, name: str, low=None, high=None, integer=False) -> pd.Series:
    """Missing tokens -> NA; any other non-numeric or out-of-range value stops the build."""
    out = []
    for value in series:
        if _is_missing(value):
            out.append(np.nan)
            continue
        try:
            number = float(value)
        except (TypeError, ValueError):
            raise ValueError(f"Non-numeric value in {name}: {value!r}") from None
        if not np.isfinite(number) or (low is not None and number < low) or (high is not None and number > high) \
                or (integer and number != int(number)):
            raise ValueError(f"Invalid value in {name}: {value!r}")
        out.append(number)
    return pd.Series(out, index=series.index, dtype="float64")


def decode_endpoint(status, raw_time, codes):
    """Return (event, months, eligibility_reason). Missing stays missing; unknown/invalid values fail."""
    missing_status, missing_time = _is_missing(status), _is_missing(raw_time)
    if not missing_status and status not in codes:
        raise ValueError(f"Unknown endpoint status: {status!r}")
    event = None if missing_status else codes[status]
    months = None if missing_time else float(raw_time)
    if months is not None and (not np.isfinite(months) or months < 0):
        raise ValueError(f"Invalid endpoint time: {raw_time!r}")
    reason = ("missing_status" if missing_status else "missing_time" if missing_time
              else "zero_time_requires_review" if months == 0 else "eligible")
    return event, months, reason


def add_endpoints(clinical: pd.DataFrame) -> pd.DataFrame:
    for prefix, (status_col, time_col, codes) in ENDPOINTS.items():
        decoded = [decode_endpoint(st, t, codes) for st, t in zip(clinical[status_col], clinical[time_col], strict=True)]
        clinical[f"{prefix}_event"] = pd.array([d[0] for d in decoded], dtype="Int8")
        clinical[f"{prefix}_months"] = pd.array([d[1] for d in decoded], dtype="Float64")
        clinical[f"{prefix}_eligibility_reason"] = [d[2] for d in decoded]
        clinical[f"{prefix}_eligible"] = clinical[f"{prefix}_eligibility_reason"].eq("eligible")
        clinical[f"{prefix}_time_unit"] = "months"
    return clinical


def build_clinical(meta: pd.DataFrame):
    missing = [c for c in REQUIRED_CHARS if c not in meta]
    if missing:
        raise ValueError(f"Missing GEO characteristics: {missing}")
    for name in ("extreme", "status.all"):  # source labels may be missing, never silently mapped
        meta[name] = meta[name].map(lambda v: v if _is_missing(v) else str(v).strip().lower())
    fusion = parse_number(meta["fusion"], "fusion", 0, 1, integer=True)
    c = pd.DataFrame({
        "sample_id": meta.geo_accession, "geo_accession": meta.geo_accession, "patient_id": meta.title,
        "platform_id": meta.platform_id,
        "age": parse_number(meta["age"], "age", 0, 120, integer=True),
        "gleason_score": parse_number(meta["gleason"], "gleason", 2, 10, integer=True),
        "gleason_primary": parse_number(meta["major.gleason"], "major.gleason", 1, 9, integer=True),
        "gleason_secondary": parse_number(meta["minor.gleason"], "minor.gleason", 1, 9, integer=True),
        "diagnosis_year": parse_number(meta["diag.yr"], "diag.yr", 1900, 2010, integer=True),
        "cancer_percent": parse_number(meta["cancer.percent"], "cancer.percent", 0, 100),
        "batch": meta["batch"], "tmprss2_erg_fusion": fusion,
        "outcome_group": meta["extreme"], "vital_status": meta["status.all"],
        "followup_months": parse_number(meta["fup.month"], "fup.month", 0),
    })
    known = c.dropna(subset=["gleason_score", "gleason_primary", "gleason_secondary"])
    if (known.gleason_primary + known.gleason_secondary != known.gleason_score).any():  # stronger than a range check
        raise ValueError("Gleason primary + secondary differs from the deposited Gleason score")
    c = add_endpoints(c)
    c["split"] = np.where(np.arange(len(c)) < N_TRAIN, "train", "validation")  # relies on deposited sample order
    c["stage"] = "T1-T2 N0 (cohort-level; per-patient stage not in GEO)"
    c["treatment"] = "none (TURP at diagnosis; watchful waiting)"
    c["tissue"] = "FFPE tumor, primary"
    return c.reset_index(drop=True)


def build_audit(meta: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    """Every GEO profile with its disposition. All 281 are single FFPE primary tumors, so all are selected."""
    return pd.DataFrame({"sample_id": clinical.sample_id, "geo_accession": meta.geo_accession,
                         "patient_id": meta.title, "platform_id": meta.platform_id, "selected": True,
                         "selection_reason": "primary_tumor",
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
    endpoint_counts = {}
    for prefix in ENDPOINTS:  # one definition of events/censored: eligible records only
        event, months, eligible = clinical[f"{prefix}_event"], clinical[f"{prefix}_months"], clinical[f"{prefix}_eligible"]
        if clinical.loc[eligible, [f"{prefix}_event", f"{prefix}_months"]].isna().any().any() or not months[eligible].gt(0).all():
            raise ValueError(f"Eligible {prefix} records need an event and a positive time")
        endpoint_counts[prefix] = {
            "eligible": int(eligible.sum()), "events": int(event[eligible].eq(1).sum()),
            "censored": int(event[eligible].eq(0).sum()),
            "median_observed_months": float(months[eligible].median()) if eligible.any() else None}
    p = PRIMARY_ENDPOINT
    counts = {"samples": len(clinical), "patients": int(clinical.patient_id.nunique()),
              "source_features": len(mapping), "retained_feature_rows": len(expression),
              "gene_rows": int(genes.gene_id.nunique()), "unique_gene_symbols": int(genes.gene_symbol.nunique()),
              "features_sharing_gene_id": int((genes.n_probes_same_gene > 1).sum()),
              f"{p}_events": endpoint_counts[p]["events"], f"{p}_censored": endpoint_counts[p]["censored"],
              "train": int((clinical.split == "train").sum()), "validation": int((clinical.split == "validation").sum())}
    odd = clinical.loc[(clinical.gleason_primary > 5) | (clinical.gleason_secondary > 5), "sample_id"].tolist()
    return {"gleason_pattern_above_5_samples": odd, "structurally_valid": True, "counts": counts, "primary_endpoint": p, "endpoints": endpoint_counts,
            "excluded_features_by_status": mapping.loc[~mapping.selected, "mapping_status"].value_counts().to_dict(),
            "expression_range": [float(expression.to_numpy().min()), float(expression.to_numpy().max())],
            "missing_clinical_values": {k: int(v) for k, v in clinical.isna().sum().items() if v},
            "sample_order_aligned": True, "feature_order_aligned": True}


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
        "endpoints": {prefix: {"status": st, "months": t, "status_codes": codes}
                      for prefix, (st, t, codes) in ENDPOINTS.items()},
        "primary_endpoint": PRIMARY_ENDPOINT,
        "endpoint_policy": {"unit": "months", "zero_time_policy": "retain but exclude from eligibility pending review",
                            "missing_policy": "retain nullable event/time and eligibility reason; never assume censored"},
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
                        "Extreme-groups design",
                        "One sample (GSM416492) has minor Gleason pattern 6 (3+6=9), outside the usual 1-5 pattern range but consistent with its deposited total; retained as deposited and listed in validation.json"],
        "code_sha256": code_hashes(__name__),
    }
    metadata_out["output_schemas"] = {name: {c: {"dtype": str(f[c].dtype), "nullable": bool(f[c].isna().any())} for c in f.columns}
                                      for name, f in [("clinical.csv", clinical), ("samples.csv", samples), ("genes.csv", genes),
                                                      ("feature_mapping.csv.gz", mapping), ("sample_audit.csv", audit)]}
    # Rebuild with the published provenance when verifying an existing output pair; every staged
    # file must still match its published counterpart byte for byte (same rule as the Taylor pipeline).
    published_metadata = output_root / snapshot / "metadata.json"
    if published_metadata.is_file() and (anndata_root / snapshot).is_dir():
        try:
            published = json.loads(published_metadata.read_text())
        except ValueError:
            published = None  # the destination check reports the mismatch
        if isinstance(published, dict) and isinstance(published.get("code_sha256"), dict):
            metadata_out["code_sha256"] = published["code_sha256"]
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
