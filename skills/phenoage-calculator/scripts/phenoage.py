#!/usr/bin/env python3
"""Blood PhenoAge — Levine 2018 biological-age clock (local, zero-dependency).

MIT License. Copyright (c) 2026 Insilico-org (Insilico Medicine) — original
verified implementation in Insilico-org/plsi, src/plsi/map/derive.py.
Portions derived and adapted under the same MIT terms. Formula and public-domain
science after Levine 2018 (doi:10.18632/aging.101414, PMC5940111) and the BioAge R
reference implementation (Kwon & Belsky 2021, doi:10.1007/s11357-021-00480-5).

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in the
Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Not biolearn's DNAm 'PhenoAge' (that is a methylation clock). This is the blood-only
PhenoAge, computed from 9 standard clinical markers + chronological age. Every value
is unit-normalized to the canonical units the formula expects.

Inputs (canonical units the formula expects):
  albumin        g/L
  creatinine     umol/L
  glucose        mmol/L
  crp            mg/dL          (ln(crp); must be > 0)
  lymphocyte_pct % (of WBC)
  mcv            fL
  rdw            % (red cell distribution width)
  alp            U/L
  wbc            10^9/L

Run read-only; nothing leaves this machine.

Usage:
  python3 phenoage.py --age 52 --json markers.json
  python3 phenoage.py --age 52 '{"glucose_mg":"85","crp_mgl":"1.1","creatinine_mgd":"0.96"}'
"""

import math
import sys
import json

# Levine 2018 coefficients (github.com/dayoonkwon/BioAge R, orig=TRUE).
# The 0.00188 ALP weight matches the furniture in the PLSI extraction; Levine's own
# table prints 0.0019 (rounding) — the BioAge R reference carries 0.00188.
WEIGHTS = {
    "albumin": -0.0336,
    "creatinine": 0.0095,
    "glucose": 0.1953,
    "lymphocyte_pct": -0.0120,
    "mcv": 0.0268,
    "rdw": 0.3306,
    "alp": 0.00188,
    "wbc": 0.0554,
}
BASE = -19.907
AGE_W = 0.0804
GAMMA = 0.0076927


def phenoage(inputs: dict, age: float) -> float:
    """Blood PhenoAge from canonical-unit inputs keyed by the names above."""
    if inputs.get("crp", 0) <= 0:
        raise ValueError("PhenoAge needs CRP > 0 (the formula takes ln CRP)")
    xb = BASE + 0.0954 * math.log(inputs["crp"]) + AGE_W * age
    xb += sum(w * inputs[k] for k, w in WEIGHTS.items())
    mortality = 1 - math.exp(-math.exp(xb) * (math.exp(120 * GAMMA) - 1) / GAMMA)
    return 141.50225 + math.log(-0.00553 * math.log(1 - mortality)) / 0.090165


# ---- unit input helpers: accept common clinical units, normalise ----
# keys accepted -> (canonical key, factor, unit name)
_GLU = 1.0 / 18.016        # mg/dL -> mmol/L
_CRP = 1.0 / 10.0          # mg/L  -> mg/dL
_CRE = 88.4                # mg/dL -> umol/L


def consume(v, raw):
    """Return (canonical_value, None) or raise with hint if key unknown."""
    k = raw["key"]
    key = raw.get("canonical") or k
    x = float(v)
    if k == "glucose_mg": return glucose_mmol(x), _GLU
    if k == "glucose_mmol": return x, None
    if k == "crp_mgl": return crp_mgd(x), _CRP
    if k == "crp_mgd": return x, None
    if k == "creatinine_mgd": return creat_umol(x), _CRE
    if k == "creatinine_umol": return x, None
    return x, None


def glucose_mmol(mg_dl): return mg_dl / 18.016
def crp_mgd(mg_l): return mg_l / 10.0
def creat_umol(mg_dl): return mg_dl * 88.4


def _resolve_unit(key):
    """Return (canonical_key, converter). Unknown/plain keys pass through."""
    if key == "glucose_mg":    return "glucose", glucose_mmol
    if key == "glucose_mmol":  return "glucose", float
    if key == "crp_mgl":       return "crp", crp_mgd
    if key == "crp_mgd":       return "crp", float
    if key == "creatinine_mgd":return "creatinine", creat_umol
    if key == "creatinine_umol": return "creatinine", float
    return key.replace("canonical:", ""), float


def main():
    args = sys.argv[1:]
    age = None
    data = None

    for x in args:
        if x.startswith("{"):
            data = json.loads(x)
        elif x.startswith("--age="):
            age = float(x.split("=", 1)[1])
    if data is None:
        try:
            data = json.load(open(args[-1]))
        except Exception:
            data = None
    if data is None:
        print("no input"); return 2

    # resolve each marker to a canonical-name float, tolerating blanks
    canonical = {}
    for raw in data.get("markers", []):
        v = raw.get("value")
        if v is None:
            continue
        ckey, conv = _resolve_unit(raw["key"])
        canonical[ckey] = conv(float(v))

    if age is None:
        d = data.get("age")
        age = float(d) if d is not None else None

    CRIT = set(WEIGHTS) | {"crp"}
    if age is not None and not (set(canonical) >= CRIT):
        out = {
            "complete": False,
            "chronological_age": age,
            "phenoage": None,
            "phenoage_acceleration": None,
            "have_markers": sorted(set(canonical) & CRIT),
            "missing_markers": sorted(CRIT - set(canonical)),
            "note": "Cannot compute a valid PhenoAge until all 9 markers are present.",
        }
        print(json.dumps(out, indent=2))
        return 0

    if age is None:
        print("no age given (pass --age or an 'age' field)"); return 2

    val = phenoage(canonical, age)
    print(json.dumps({
        "complete": True,
        "chronological_age": age,
        "phenoage": round(val, 2),
        "phenoage_acceleration": round(val - age, 2),
        "markers_used": {k: canonical[k] for k in sorted(CRIT)},
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
