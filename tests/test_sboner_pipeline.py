import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

from ds_bu_bioinformatics import sboner_pipeline as sp

N = 4
GSM = [f"GSM{i}" for i in range(1, N + 1)]


def matrix_text(extreme=("Lethal", "Indolent", "Lethal", "Indolent"), status=("Dead", "Alive", "Dead", "Dead"),
                fusion=("0", "1", "NA", "0"), probes=("P1", "P2", "P3")):
    q = lambda vals: "\t".join(f'"{v}"' for v in vals)
    lines = [f'!Series_geo_accession\t"{sp.ACC}"',
             f"!Sample_title\t{q([f'prostate_{i}' for i in range(N)])}",
             f"!Sample_geo_accession\t{q(GSM)}",
             f"!Sample_platform_id\t{q([sp.PLATFORM] * N)}",
             f"!Sample_data_processing\t{q(['log2 qspline'] * N)}"]
    chars = {"gleason": ["7", "6", "9", "7"], "major.gleason": ["4", "3", "4", "3"], "minor.gleason": ["3", "3", "5", "4"],
             "batch": ["1", "1", "2", "2"], "extreme": list(extreme), "cancer.percent": ["20", "30", "NA", "15"],
             "age": ["61", "80", "86", "74"], "diag.yr": ["1996", "1991", "1995", "1988"],
             "status.all": list(status), "fup.month": ["122", "124", "54", "149"], "fusion": list(fusion)}
    for k, v in chars.items():
        lines.append("!Sample_characteristics_ch1\t" + q([f"{k}: {x}" for x in v]))
    lines += ["!series_matrix_table_begin", "ID_REF\t" + q(GSM)]
    for i, p in enumerate(probes):
        lines.append(f'"{p}"\t' + "\t".join(f"{8 + i + j * 0.123456789012345:.15g}" for j in range(N)))
    lines.append("!series_matrix_table_end")
    return "\n".join(lines) + "\n"


def gpl_text(rows):
    head = f"^PLATFORM = {sp.PLATFORM}\n!platform_table_begin\nID\tGeneID\tREFSEQ\tGB_ACC\tSymbol\tDescription\n"
    return head + "".join("\t".join(r) + "\n" for r in rows) + "!platform_table_end\n"


GOOD_GPL = [("P1", "100", "NM_000001 ; NP_000001", "NM_000001", "AAA", "a"),
            ("P2", "100", "NM_000002 ; NP_000002", "NM_000002", "AAA", "a2"),  # shares gene with P1
            ("P3", "300", "NM_000003", "NM_000003", "CCC", "c")]


def write_gz(path: Path, text: str):
    with gzip.open(path, "wt") as f:
        f.write(text)


class StubCache:
    def __init__(self, root: Path, matrix: str, gpl: str):
        self.root, self.manifest_path = root, root / "manifest.json"
        write_gz(root / sp.MATRIX_NAME, matrix)
        write_gz(root / sp.GPL_NAME, gpl)
        self.manifest_path.write_text("{}")

    def verify_all(self):
        pass

    def fetch(self, name, url):
        return self.root / name


class SbonerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def parse(self, matrix):
        p = self.root / "m.gz"
        write_gz(p, matrix)
        return sp.read_geo(p)

    def test_read_geo_preserves_order_and_values(self):
        meta, expr, _ = self.parse(matrix_text())
        self.assertEqual(list(expr.columns), GSM)
        self.assertEqual(list(expr.index), ["P1", "P2", "P3"])
        self.assertEqual(meta.loc[0, "gleason"], "7")

    def test_read_geo_rejects_bad_inputs(self):
        with self.assertRaises(ValueError):
            self.parse(matrix_text().replace(sp.ACC, "GSE1"))
        with self.assertRaises(ValueError):  # duplicate probe
            self.parse(matrix_text(probes=("P1", "P1", "P3")))
        with self.assertRaises(ValueError):  # truncated table
            self.parse(matrix_text().replace("!series_matrix_table_end\n", ""))
        with self.assertRaises(ValueError):  # non-finite expression
            self.parse(matrix_text().replace("8.", "nan.", 1).replace("nan.", "nan"))

    def test_mapping_keeps_features_sharing_a_gene(self):
        p = self.root / "g.gz"
        write_gz(p, gpl_text(GOOD_GPL))
        m = sp.resolve_mapping(sp.parse_platform(p), ["P1", "P2", "P3"])
        self.assertTrue(m.selected.all())
        expr = pd.DataFrame(np.arange(12.0).reshape(3, 4), index=pd.Index(["P1", "P2", "P3"]), columns=GSM)
        out, genes = sp.select_features(expr, m)
        self.assertEqual(len(out), 3)  # no averaging
        self.assertEqual(genes.gene_id.tolist(), ["100", "100", "300"])
        self.assertEqual(genes.n_probes_same_gene.tolist(), [2, 2, 1])

    def test_mapping_excludes_invalid_and_conflicting(self):
        bad = GOOD_GPL + [("P4", "", "NM_000004", "NM_000004", "DDD", ""),
                          ("P5", "500", "NM_000005", "NM_000005", "", ""),
                          ("P6", "600", "XX_1", "XX_1", "FFF", ""),
                          ("P7", "700", "NM_000007", "NM_000007", "AAA", "")]  # symbol AAA on two gene IDs
        p = self.root / "g.gz"
        write_gz(p, gpl_text(bad))
        m = sp.resolve_mapping(sp.parse_platform(p), [r[0] for r in bad]).set_index("feature_id")
        self.assertEqual(m.loc["P4", "mapping_status"], "no_entrez_id")
        self.assertEqual(m.loc["P5", "mapping_status"], "no_symbol")
        self.assertEqual(m.loc["P6", "mapping_status"], "invalid_refseq")
        self.assertEqual(m.loc["P7", "mapping_status"], "gene_symbol_identity_conflict")
        self.assertEqual(m.loc["P1", "mapping_status"], "gene_symbol_identity_conflict")

    def test_mapping_rejects_feature_set_mismatch(self):
        p = self.root / "g.gz"
        write_gz(p, gpl_text(GOOD_GPL))
        with self.assertRaises(ValueError):
            sp.resolve_mapping(sp.parse_platform(p), ["P1", "P2"])

    def test_clinical_events_and_nullables(self):
        meta, _, _ = self.parse(matrix_text())
        c = sp.build_clinical(meta)
        self.assertEqual(c.pcsm_event.tolist(), [1, 0, 1, 0])
        self.assertEqual(c.os_event.tolist(), [1, 0, 1, 1])
        for prefix in ("pcsm", "os"):  # shared endpoint contract
            for suffix in ("event", "months", "eligible", "eligibility_reason", "time_unit"):
                self.assertIn(f"{prefix}_{suffix}", c)
        self.assertTrue(c.pcsm_eligible.all())
        self.assertEqual(c.pcsm_months.tolist(), [122, 124, 54, 149])
        self.assertEqual(set(c.pcsm_time_unit), {"months"})
        self.assertTrue(np.isnan(c.tmprss2_erg_fusion[2]))  # NA stays missing, not 0
        self.assertTrue(np.isnan(c.cancer_percent[2]))
        self.assertEqual(c.split.tolist(), ["train"] * 4)

    def test_missing_outcome_stays_missing_and_ineligible(self):
        meta, _, _ = self.parse(matrix_text(extreme=("Lethal", "NA", "Lethal", "Indolent")))
        meta.loc[3, "fup.month"] = "0"
        c = sp.build_clinical(meta)
        self.assertTrue(pd.isna(c.pcsm_event[1]))  # never assumed censored
        self.assertEqual(c.pcsm_eligibility_reason.tolist(),
                         ["eligible", "missing_status", "eligible", "zero_time_requires_review"])
        self.assertEqual(c.pcsm_eligible.tolist(), [True, False, True, False])

    def test_strict_covariates_fail_on_bad_values(self):
        for column, bad in [("age", "200"), ("gleason", "11"), ("gleason", "6.5"), ("fup.month", "-3"),
                            ("fusion", "2"), ("cancer.percent", "abc")]:
            meta, _, _ = self.parse(matrix_text())
            meta.loc[0, column] = bad
            with self.assertRaises(ValueError, msg=f"{column}={bad}"):
                sp.build_clinical(meta)

    def test_audit_uses_shared_columns(self):
        meta, _, _ = self.parse(matrix_text())
        audit = sp.build_audit(meta, sp.build_clinical(meta))
        self.assertTrue({"sample_id", "selected", "selection_reason"} <= set(audit.columns))

    def test_clinical_rejects_unknown_codes(self):
        meta, _, _ = self.parse(matrix_text(extreme=("Lethal", "Maybe", "Lethal", "Indolent")))
        with self.assertRaises(ValueError):
            sp.build_clinical(meta)

    def test_end_to_end_build_and_protection(self):
        raw = self.root / "raw"
        raw.mkdir()
        cache = StubCache(raw, matrix_text(), gpl_text(GOOD_GPL))
        roots = [self.root / n for n in ("processed", "interim", "anndata")]
        with mock.patch.object(sp, "EXPECTED_SAMPLES", N), mock.patch.object(sp, "N_TRAIN", 3):
            checks = sp.build(cache, *roots, snapshot="t")
            self.assertEqual(checks["counts"]["retained_feature_rows"], 3)
            self.assertTrue(checks["structurally_valid"])
            self.assertEqual(checks["primary_endpoint"], "pcsm")
            self.assertEqual(checks["counts"]["pcsm_events"], 2)
            self.assertEqual(checks["endpoints"]["pcsm"]["eligible"], 4)
            self.assertEqual(checks["counts"]["validation"], 1)
            self.assertTrue(checks["anndata_round_trip_valid"])
            sp.build(cache, *roots, snapshot="t")  # identical rebuild accepted
            (roots[0] / "t" / "clinical.csv").write_text("tampered")
            with self.assertRaises(ValueError):  # differing outputs protected
                sp.build(cache, *roots, snapshot="t")
        meta = json.loads((roots[0] / "t" / "metadata.json").read_text())
        self.assertEqual(meta["geo_series"], sp.ACC)


if __name__ == "__main__":
    unittest.main()
