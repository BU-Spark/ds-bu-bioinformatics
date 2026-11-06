"""Offline tests on small synthetic data shaped like the real cBioPortal files.

Run from the repo root:  pytest -q
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import tcga_prad_pilot as p  # noqa: E402

PATIENTS = [f"TCGA-AA-{i:04d}" for i in range(1, 31)]


@pytest.fixture
def expr_raw():
    rng = np.random.default_rng(0)
    samples = [pid + "-01" for pid in PATIENTS] + ["TCGA-AA-0001-11", "TCGA-AA-0002-01B"]
    genes = ["AR", "KLK3", "ERG", "PTEN", "AR", "?", "GAPGENE"] + [f"G{i}" for i in range(40)]
    df = pd.DataFrame(rng.gamma(2, 200, (len(genes), len(samples))), columns=samples)
    df.loc[genes.index("GAPGENE"), samples[:3]] = np.nan          # a batch-correction gap
    df.insert(0, "Entrez_Gene_Id", range(len(genes)))
    df.insert(0, "Hugo_Symbol", genes)
    return df


@pytest.fixture
def patients():
    rows = {}
    for i, pid in enumerate(PATIENTS):
        rows[pid] = {"PFS_STATUS": f"{int(i % 3 == 0)}:X", "PFS_MONTHS": str(5 + i),
                     "OS_STATUS": "0:LIVING", "OS_MONTHS": str(30 + i), "AGE": str(50 + i),
                     "PATH_T_STAGE": "T2C", "PATH_N_STAGE": "N0"}
    return pd.DataFrame.from_dict(rows, orient="index")


def test_clean_expression(expr_raw):
    log = p.QCLog(verbose=False)
    expr, features, removed = p.clean_expression(expr_raw, log)
    assert "?" not in expr.index                       # bad symbol dropped
    assert expr.index.is_unique                        # duplicate AR collapsed
    assert removed == ["GAPGENE"]                      # gene with missing values removed
    assert not expr.isna().any().any()
    assert list(features.index) == list(expr.index)


def test_select_samples_barcode_rule(expr_raw):
    log = p.QCLog(verbose=False)
    expr, _, _ = p.clean_expression(expr_raw, log)
    man = p.select_samples(expr, log)
    reasons = man.set_index("sample_id")["reason"]
    assert reasons["TCGA-AA-0001-11"].startswith("not primary")
    assert reasons["TCGA-AA-0002-01B"] == "extra sample for same patient"
    assert man.loc[man["included"], "patient_id"].is_unique


def test_match_and_outputs(expr_raw, patients, tmp_path):
    log = p.QCLog(verbose=False)
    expr, features, removed = p.clean_expression(expr_raw, log)
    man = p.select_samples(expr, log)
    clin = p.clean_clinical(patients, pd.DataFrame(), None, log)
    clinical, expression, man = p.match(man, clin, expr, log)
    assert list(expression.columns) == list(clinical.index)
    assert clinical["patient_id"].is_unique
    assert set(clinical["PFI_event"].unique()) <= {0, 1}
    features = p.annotate_features(features, expression)
    p.save_outputs(tmp_path, expression, clinical, features, man, {"sources": {}}, removed, log)
    for name in ["expression_log2_genes_x_samples.csv.gz", "clinical_by_sample.csv",
                 "feature_annotation.csv", "sample_manifest.csv", "data_dictionary.csv",
                 "genes_removed_missing_values.csv", "provenance.json", "qc_summary.md"]:
        assert (tmp_path / name).exists(), name


def test_module_coverage(expr_raw, tmp_path):
    log = p.QCLog(verbose=False)
    expr, _, removed = p.clean_expression(expr_raw, log)
    mod_file = tmp_path / "modules.csv"
    pd.DataFrame({"module": ["M1", "M1", "M2"], "gene": ["AR", "GAPGENE", "NOTAGENE"]}).to_csv(
        mod_file, index=False)
    cov = p.module_coverage(str(mod_file), expr.index, removed, log)
    assert cov.loc["M1", "found"] == 1 and cov.loc["M1", "lost_to_missing_values"] == 1
    assert cov.loc["M2", "found"] == 0
