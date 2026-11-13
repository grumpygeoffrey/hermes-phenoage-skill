#!/usr/bin/env python3
"""Tests for phenoage.py — stdlib only: python3 -m unittest discover -s tests"""

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from datetime import date

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "scripts", "phenoage.py")
if not os.path.exists(SCRIPT):                      # repo layout: script in skills/…/scripts
    SCRIPT = os.path.join(ROOT, "skills", "phenoage-calculator", "scripts", "phenoage.py")

spec = importlib.util.spec_from_file_location("phenoage", SCRIPT)
pa = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pa)

PANEL = {"albumin": 44, "creatinine": 88, "glucose": 5.1, "crp": 0.16,
         "lymphocyte_pct": 28.5, "mcv": 89, "rdw": 13.8, "alp": 72, "wbc": 6.2}


class TestClock(unittest.TestCase):
    def test_complete_panel_computes(self):
        v = pa.phenoage(PANEL, 52)
        self.assertTrue(40 < v < 60)

    def test_monotonic_in_crp(self):
        low = pa.phenoage({**PANEL, "crp": 0.05}, 52)
        high = pa.phenoage({**PANEL, "crp": 0.5}, 52)
        self.assertLess(low, high)

    def test_crp_zero_rejected(self):
        with self.assertRaises(ValueError):
            pa.phenoage({**PANEL, "crp": 0}, 52)

    def test_missing_marker_rejected(self):
        incomplete = dict(PANEL)
        incomplete.pop("rdw")
        with self.assertRaises(ValueError):
            pa.phenoage(incomplete, 52)


class TestUnits(unittest.TestCase):
    def test_glucose_mg_dl(self):
        self.assertAlmostEqual(pa.to_canonical("glucose", 92, "mg/dL"), 92 / 18.016, 6)

    def test_albumin_g_dl(self):
        self.assertAlmostEqual(pa.to_canonical("albumin", 4.4, "g/dL"), 44.0, 9)

    def test_creatinine_mg_dl(self):
        self.assertAlmostEqual(pa.to_canonical("creatinine", 1.0, "mg/dL"), 88.4, 9)

    def test_crp_mg_l(self):
        self.assertAlmostEqual(pa.to_canonical("crp", 1.6, "mg/L"), 0.16, 9)

    def test_wbc_per_ul(self):
        self.assertAlmostEqual(pa.to_canonical("wbc", 6200, "/uL"), 6.2, 9)

    def test_wbc_10_3_ul(self):
        self.assertAlmostEqual(pa.to_canonical("wbc", 6.2, "10^3/uL"), 6.2, 9)

    def test_alp_ukat(self):
        self.assertAlmostEqual(pa.to_canonical("alp", 1.2, "ukat/L"), 72.0, 9)

    def test_generic_ucum_g_l(self):
        # mg/dL -> g/L via prefix/volume scaling
        self.assertAlmostEqual(pa.to_canonical("albumin", 440, "mg/dL"), 4.4, 9)

    def test_incompatible_unit_rejected(self):
        with self.assertRaises(ValueError):
            pa.to_canonical("glucose", 5, "kg")

    def test_censored_rejected(self):
        with self.assertRaises(ValueError) as ctx:
            pa.to_canonical("crp", "<0.03", "mg/dL")
        self.assertIn("censored", str(ctx.exception))


class TestAge(unittest.TestCase):
    def test_exact_age(self):
        self.assertAlmostEqual(pa.age_at_date(date(1974, 6, 15), date(2026, 11, 13)),
                               52.41, 2)

    def test_year_approximation_within_half_year(self):
        exact = pa.age_at_date(date(1974, 6, 15), date(2026, 11, 13))
        approx = pa.chronological_age(1974, date(2026, 11, 13))
        self.assertLess(abs(exact - approx), 0.5)

    def test_leap_day_birth(self):
        self.assertAlmostEqual(pa.age_at_date(date(2000, 2, 29), date(2026, 3, 1)), 26.0, 1)


class TestDrivers(unittest.TestCase):
    def test_effects_sum_to_total(self):
        later = {**PANEL, "crp": 0.11, "rdw": 13.2, "albumin": 45}
        ch = pa.explain_change({"values": PANEL, "age": 51.4, "id": "a"},
                               {"values": later, "age": 52.4, "id": "b"})
        self.assertAlmostEqual(sum(e["years"] for e in ch["effects"]), ch["total"], 12)
        self.assertAlmostEqual(ch["total"],
                               pa.phenoage(later, 52.4) - pa.phenoage(PANEL, 51.4), 12)

    def test_identical_panels_only_age(self):
        ch = pa.explain_change({"values": PANEL, "age": 51.4},
                               {"values": dict(PANEL), "age": 51.4})
        self.assertAlmostEqual(ch["total"], 0.0, 12)
        self.assertFalse([e for e in ch["effects"] if e["kind"] == "marker"])

    def test_marker_named_with_direction(self):
        worse = {**PANEL, "crp": 0.40}
        ch = pa.explain_change({"values": PANEL, "age": 52}, {"values": worse, "age": 52})
        effect = next(e for e in ch["effects"] if e.get("marker") == "crp")
        self.assertGreater(effect["years"], 0)
        self.assertAlmostEqual(ch["total"], effect["years"], 12)


class TestCli(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, SCRIPT, *args],
                              capture_output=True, text=True)

    def test_json_panel(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump({"age": 52, "markers": PANEL}, fh)
            path = fh.name
        try:
            r = self.run_cli("--format", "json", path)
            self.assertEqual(r.returncode, 0, r.stderr)
            out = json.loads(r.stdout)
            self.assertTrue(out["complete"])
            self.assertAlmostEqual(out["phenoage_acceleration"],
                                   out["phenoage"] - 52, 1)
        finally:
            os.unlink(path)

    def test_history_csv_with_mixed_units(self):
        csv = ("date,birth_year,albumin_g_dl,creatinine_mg_dl,glucose_mg_dl,"
               "crp_mg_l,lymphocyte_pct,mcv_fl,rdw_pct,alp_u_l,wbc_10_3_ul\n"
               "2025-11-20,1974,4.4,1.00,92,1.6,28.5,89.0,13.8,72,6.2\n"
               "2026-11-13,1974,4.5,0.96,90,1.1,30.1,88.5,13.2,68,5.9\n")
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as fh:
            fh.write(csv)
            path = fh.name
        try:
            r = self.run_cli("--history", path, "--format", "json")
            self.assertEqual(r.returncode, 0, r.stderr)
            out = json.loads(r.stdout)
            self.assertEqual(len(out["draws"]), 2)
            self.assertLess(out["changes"][0]["total"], 0)
        finally:
            os.unlink(path)

    def test_missing_marker_is_reported_not_invented(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump({"age": 52, "markers": {"glucose": 5.1, "crp": 0.16}}, fh)
            path = fh.name
        try:
            r = self.run_cli("--format", "json", path)
            self.assertEqual(r.returncode, 0)
            out = json.loads(r.stdout)
            self.assertFalse(out["complete"])
            self.assertIsNone(out["phenoage"])
            self.assertIn("albumin", out["missing_markers"])
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
