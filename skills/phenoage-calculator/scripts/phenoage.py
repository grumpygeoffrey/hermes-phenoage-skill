#!/usr/bin/env python3
"""
Blood PhenoAge — Levine 2018 biological-age clock (local, zero-dependency).

MIT License. Copyright (c) 2026 Insilico-org (Insilico Medicine) — original
verified implementation in Insilico-org/plsi, src/plsi/map/derive.py.
Driver decomposition and UCUM unit scaling adapted from
brainarts-lab/hebe-longevity (MIT, Copyright (c) 2026 BrainArts), which in turn
vendors PLSI under the same MIT terms. Formula and public-domain science after
Levine 2018 (doi:10.18632/aging.101414, PMC5940111) and the BioAge R reference
implementation (Kwon & Belsky 2021, doi:10.1007/s11357-021-00480-5).

Not biolearn's DNAm 'PhenoAge' (that is a methylation clock). This is the
blood-only PhenoAge, computed from 9 standard clinical markers + age.

Canonical units the formula expects:
  albumin        g/L
  creatinine     umol/L
  glucose        mmol/L
  crp            mg/dL          (ln(crp); must be > 0)
  lymphocyte_pct % (of WBC)
  mcv            fL
  rdw            % (red cell distribution width)
  alp            U/L
  wbc            10^9/L

What is new in v2.0
  1. explain_change()  — which marker moved PhenoAge between two draws, in years.
                         Single-marker substitution; a 'combined' residual absorbs
                         the non-linearity so the effects always sum to the total.
  2. chronological_age()/age_at_date() — exact age at the draw from birth year/date.
  3. --history draws.csv — trend across many draws + drivers per step.
  4. UCUM-style unit scaling — mg/dL -> g/L, 10^3/uL -> 10^9/L, % etc., plus
                         explicit conversions (glucose, creatinine, CRP, ALP).

Run read-only; nothing leaves this machine.

Usage:
  python3 phenoage.py --age 52 markers.json
  python3 phenoage.py --birth-year 1974 --json '{"markers": {"crp": 1.1, ...}}'
  python3 phenoage.py --history draws.csv --birth-year 1974
  python3 phenoage.py --history draws.csv --birth-date 1974-06-15 --format json
"""

import calendar
import csv
import json
import math
import re
import sys
from datetime import date, datetime

# ---------------------------------------------------------------- the clock

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
CRIT = set(WEIGHTS) | {"crp"}

CANON_UNITS = {
    "albumin": "g/L",
    "creatinine": "umol/L",
    "glucose": "mmol/L",
    "crp": "mg/dL",
    "lymphocyte_pct": "%",
    "mcv": "fL",
    "rdw": "%",
    "alp": "U/L",
    "wbc": "10^9/L",
}


def phenoage(inputs, age):
    """Blood PhenoAge from canonical-unit inputs keyed by the names above."""
    if inputs.get("crp", 0) <= 0:
        raise ValueError(
            "PhenoAge needs CRP > 0 (the formula takes ln CRP). "
            "A lab report of '<0.3' cannot be used as a number."
        )
    missing = CRIT - set(inputs)
    if missing:
        raise ValueError("missing markers: " + ", ".join(sorted(missing)))
    xb = BASE + 0.0954 * math.log(inputs["crp"]) + AGE_W * age
    xb += sum(w * inputs[k] for k, w in WEIGHTS.items())
    mortality = 1 - math.exp(-math.exp(xb) * (math.exp(120 * GAMMA) - 1) / GAMMA)
    return 141.50225 + math.log(-0.00553 * math.log(1 - mortality)) / 0.090165


# ------------------------------------------------------------------- age

def chronological_age(birth_year, at):
    """Age at the draw, PLSI-compatible approximation (assumes a mid-year birth).

    year_fraction = year + (day_of_year - 1) / days_in_year.
    Accurate to +/- 0.5 year; use age_at_date() when the birth date is known.
    """
    days = 366 if calendar.isleap(at.year) else 365
    return at.year + (at.timetuple().tm_yday - 1) / days - (birth_year + 0.5)


def age_at_date(birth_date, at):
    """Exact age in years (365.25-day years) at a given date."""
    years = at.year - birth_date.year
    if (at.month, at.day) < (birth_date.month, birth_date.day):
        years -= 1
    try:
        anniversary = birth_date.replace(year=birth_date.year + years)
    except ValueError:                     # 29 February
        anniversary = birth_date.replace(year=birth_date.year + years, day=28)
    return years + (at - anniversary).days / 365.25


# ---------------------------------------------------------------- units

_UNIT_EXPLICIT = {
    "albumin": {"g/l": 1.0, "g/dl": 10.0},
    "creatinine": {"umol/l": 1.0, "mg/dl": 88.4, "mg/l": 8.84},
    "glucose": {"mmol/l": 1.0, "mg/dl": 1 / 18.016, "g/l": 5.5509},
    "crp": {"mg/dl": 1.0, "mg/l": 0.1},
    "lymphocyte_pct": {"%": 1.0, "fraction": 100.0, "0-1": 100.0},
    "mcv": {"fl": 1.0, "um3": 1.0, "um^3": 1.0},
    "rdw": {"%": 1.0},
    "alp": {"u/l": 1.0, "iu/l": 1.0, "ukat/l": 60.0},
    "wbc": {
        "10^9/l": 1.0, "10*9/l": 1.0, "k/ul": 1.0, "10^3/ul": 1.0,
        "/ul": 0.001, "/mm3": 0.001, "10^3/mm3": 1.0,
    },
}

_PREFIX = {"": 1.0, "k": 1e3, "d": 1e-1, "c": 1e-2, "m": 1e-3,
           "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}
_BASES = ("[iu]", "mol", "kat", "eq", "g", "l", "u")


def _norm_unit(u):
    if u is None:
        return None
    u = (u.strip().lower()
         .replace("µ", "u").replace("μ", "u").replace("μ", "u"))
    return re.sub(r"\s+", "", u)


def _term(t):
    """One side of a UCUM unit -> (base, scale). mg -> (g, 1e-3); 10^3 -> (#, 1e3)."""
    if m := re.fullmatch(r"10[\*^](\d+)", t):
        return "#", 10.0 ** int(m.group(1))
    if t in ("", "1"):
        return "", 1.0
    if t == "%":
        return "", 0.01
    t = re.sub(r"(?<!\[)iu(?!\])", "[iu]", t)
    for base in _BASES:
        if t.endswith(base) and t[: -len(base)] in _PREFIX:
            return base, _PREFIX[t[: -len(base)]]
    return None


def _scale(unit):
    num, _, den = unit.partition("/")
    a, b = _term(num), (_term(den) if den else ("", 1.0))
    if a and b and "/" not in den:
        return (a[0], b[0]), a[1] / b[1]
    return None


def to_canonical(marker, value, unit=None, _warn=None):
    """Convert a value to the marker's canonical unit.

    Order: identity -> explicit clinical conversion -> generic UCUM scale.
    Raises ValueError rather than guessing between kinds of quantity.
    """
    if marker not in CANON_UNITS:
        raise ValueError(f"unknown marker {marker!r}")
    target = CANON_UNITS[marker]
    sv = str(value).strip()
    if sv[:1] in ("<", ">"):
        raise ValueError(
            f"{marker} reported as {sv!r}: censored values are not supported. "
            "The formula needs a number; substituting the bound would silently "
            "bias the result, so get the numeric value or leave the marker out."
        )
    if unit is None:
        return float(value)
    u = _norm_unit(unit)
    t = _norm_unit(target)
    if u == t:
        return float(value)
    factor = _UNIT_EXPLICIT.get(marker, {}).get(u)
    if factor is not None:
        return float(value) * factor
    a, b = _scale(u), _scale(t)
    if a and b and a[0] == b[0]:
        f = a[1] / b[1]
        if _warn is not None and marker in ("albumin", "crp"):
            pass
        return float(value) * f
    raise ValueError(
        f"cannot convert {marker} {unit!r} to {target!r} "
        "(no clinical conversion or matching UCUM base)"
    )


# Legacy JSON keys used by v1, plus canonical names.
_KEY_UNITS = {
    "glucose_mg": ("glucose", "mg/dL"),
    "glucose_mmol": ("glucose", "mmol/L"),
    "crp_mgl": ("crp", "mg/L"),
    "crp_mgd": ("crp", "mg/dL"),
    "creatinine_mgd": ("creatinine", "mg/dL"),
    "creatinine_umol": ("creatinine", "umol/L"),
    "albumin_gdl": ("albumin", "g/dL"),
}


def resolve_key(key):
    """Map a v1-style key ('glucose_mg') or a canonical name to (marker, unit)."""
    k = key.replace("canonical:", "").strip()
    if k in _KEY_UNITS:
        return _KEY_UNITS[k]
    if k in CANON_UNITS:
        return k, CANON_UNITS[k]
    raise ValueError(f"unknown marker key {key!r}")


# ------------------------------------------------------------------ drivers

def explain_change(previous, current):
    """What moved PhenoAge between two draws, in years.

    Method: start from the previous draw, replace one marker at a time with its
    new value, recompute. The formula is non-linear, so the single effects do not
    add up exactly; the remainder is reported as 'combined' and the effects always
    sum to the real total.

    previous/current: {'values': {9 canonical markers}, 'age': float, 'id': str}
    """
    before = dict(previous["values"])
    age_before = float(previous["age"])
    base = phenoage(before, age_before)
    total = phenoage(current["values"], float(current["age"])) - base

    effects = []
    for marker in sorted(before):
        after = current["values"].get(marker)
        if after is None or after == before[marker]:
            continue
        changed = dict(before)
        changed[marker] = after
        effects.append({
            "kind": "marker",
            "marker": marker,
            "from": before[marker],
            "to": after,
            "unit": CANON_UNITS[marker],
            "years": phenoage(changed, age_before) - base,
        })

    age_effect = phenoage(before, float(current["age"])) - base
    if abs(age_effect) > 1e-12:
        effects.append({
            "kind": "age",
            "from": round(age_before, 2),
            "to": round(float(current["age"]), 2),
            "years": age_effect,
        })

    rest = total - sum(e["years"] for e in effects)
    effects.sort(key=lambda e: -abs(e["years"]))
    effects.append({"kind": "combined", "years": rest})

    return {
        "from_id": previous.get("id"),
        "to_id": current.get("id"),
        "from_date": previous.get("measured_at"),
        "to_date": current.get("measured_at"),
        "total": total,
        "effects": effects,
    }


def previous_draw(history, current):
    """The draw immediately before `current` in (measured_at, id) order."""
    key = (current.get("measured_at") or "", str(current.get("id") or ""))
    earlier = [r for r in history
               if ((r.get("measured_at") or ""), str(r.get("id") or "")) < key]
    return earlier[-1] if earlier else None


# ------------------------------------------------------------------ input

_DATE_COLS = {"date", "measured_at", "draw_date", "datum", "afnamedatum"}
_YEAR_COLS = {"birth_year", "geboortejaar", "year_of_birth"}
_AGE_COLS = {"age", "leeftijd"}
_ID_COLS = {"id", "draw_id", "sample_id"}

_UNIT_TOKENS = {
    "g_l": "g/L", "g_dl": "g/dL",
    "umol_l": "umol/L", "mg_dl": "mg/dL", "mg_l": "mg/L",
    "mmol_l": "mmol/L", "fl": "fL", "um3": "um3", "um_3": "um3",
    "u_l": "U/L", "iu_l": "IU/L", "ukat_l": "ukat/L",
    "109_l": "10^9/L", "10_9_l": "10^9/L", "10_3_ul": "10^3/uL",
    "k_ul": "K/uL", "ul": "/uL", "mm3": "/mm3", "pct": "%",
}


def _norm_header(h):
    h = h.strip().lower().replace("µ", "u").replace("μ", "u")
    h = h.replace("%", "pct").replace("^", "").replace("/", "_")
    h = re.sub(r"[\s\-]+", "_", h)
    return re.sub(r"_+", "_", h).strip("_")


def parse_header(name):
    """CSV header -> (marker|'__date__'|..., unit). None when not a known column."""
    n = _norm_header(name)
    if n in _DATE_COLS:
        return "__date__", None
    if n in _YEAR_COLS:
        return "__birth_year__", None
    if n in _AGE_COLS:
        return "__age__", None
    if n in _ID_COLS:
        return "__id__", None
    if n in ("crp_hs", "hscrp", "crp_hoog_sensitief"):
        return "crp", "mg/L"
    for marker in sorted(CANON_UNITS, key=len, reverse=True):
        if n == marker:
            return marker, CANON_UNITS[marker]
        if n.startswith(marker + "_"):
            unit = _UNIT_TOKENS.get(n[len(marker) + 1:])
            if unit:
                return marker, unit
    return None, None


def _parse_date(value):
    v = value.strip()
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", v):
        return date.fromisoformat(v)
    if re.fullmatch(r"\d{4}-\d{2}", v):
        return date.fromisoformat(v + "-15")
    if re.fullmatch(r"\d{2}/\d{2}/\d{4}", v):
        d, m, y = v.split("/")
        return date(int(y), int(m), int(d))
    if re.fullmatch(r"\d{4}", v):
        return date(int(v), 7, 1)
    raise ValueError(f"unrecognised date {value!r}")


def _panel_from_json(data):
    """Accept {'markers': [{key,value,unit?}]}, {'markers': {name: value}}, or flat dict."""
    out = {}
    raw_markers = data.get("markers")
    if isinstance(raw_markers, dict):
        items = [{"key": k, "value": v} for k, v in raw_markers.items()]
    elif isinstance(raw_markers, list):
        items = raw_markers
    else:                                  # flat {glucose: 5.0, ...}
        items = [{"key": k, "value": v} for k, v in data.items()
                 if k in CANON_UNITS or k in _KEY_UNITS]
    for item in items:
        if item.get("value") in (None, ""):
            continue
        marker, unit = resolve_key(item["key"])
        unit = item.get("unit") or unit
        out[marker] = to_canonical(marker, str(item["value"]).strip(), unit)
    return out


def _panel_from_row(row, columns, birth_year=None, birth_date=None):
    """Build one draw from a CSV row."""
    values, at, age, ident = {}, None, None, None
    for col, (kind, unit) in columns.items():
        v = (row.get(col) or "").strip()
        if not v:
            continue
        if kind == "__date__":
            at = _parse_date(v)
        elif kind == "__birth_year__":
            birth_year = int(v)
        elif kind == "__age__":
            age = float(v)
        elif kind == "__id__":
            ident = v
        else:
            values[kind] = to_canonical(kind, v, unit)

    if at is None:
        raise ValueError("history CSV needs a date column (date or measured_at)")
    if birth_date is not None:
        age = age_at_date(birth_date, at)
    elif birth_year is not None:
        age = chronological_age(birth_year, at)
    if age is None:
        raise ValueError("no age: pass --birth-year/--birth-date or add a birth_year/age column")

    missing = CRIT - set(values)
    if missing:
        raise ValueError("draw missing markers: " + ", ".join(sorted(missing)))

    return {
        "id": ident or at.isoformat(),
        "measured_at": at.isoformat(),
        "date": at,
        "age": age,
        "values": values,
        "phenoage": phenoage(values, age),
    }


# ------------------------------------------------------------------ output

def _fmt(v, nd=2):
    return f"{v:.{nd}f}"


def print_panel(res, source=None):
    if not res["complete"]:
        print(json.dumps(res, indent=2))
        return
    print(f"Blood PhenoAge{(' — ' + source) if source else ''}")
    print(f"  chronological age : {_fmt(res['chronological_age'])} y")
    print(f"  PhenoAge          : {_fmt(res['phenoage'])} y")
    print(f"  acceleration      : {res['phenoage_acceleration']:+.2f} y")
    print("  markers (canonical units):")
    for k in sorted(res["markers_used"]):
        print(f"    {k:<15} {res['markers_used'][k]:>8.3f}  {CANON_UNITS[k]}")


def print_trend(draws, changes):
    print(f"Blood PhenoAge — {len(draws)} draws\n")
    print(f"{'DATE':<12} {'AGE':>6} {'PHENOAGE':>9} {'ACCEL':>8}")
    for d in draws:
        print(f"{d['measured_at']:<12} {d['age']:>6.1f} "
              f"{d['phenoage']:>9.2f} {d['phenoage']-d['age']:>+8.2f}")
    for ch in changes:
        print(f"\nChange {ch['from_date']} -> {ch['to_date']}: {ch['total']:+.2f} y")
        for e in ch["effects"]:
            if e["kind"] == "marker":
                print(f"  {e['marker']:<14} {e['from']:.2f} -> {e['to']:.2f} "
                      f"{e['unit']:<8} {e['years']:+.2f} y")
            elif e["kind"] == "age":
                print(f"  {'age':<14} {e['from']:.1f} -> {e['to']:.1f}   "
                      f"{'':<8} {e['years']:+.2f} y")
            else:
                print(f"  {'combined':<14} {'(non-linearity)':<24} {e['years']:+.2f} y")


# --------------------------------------------------------------------- main

def _load_json_arg(arg):
    if arg is True:
        return None
    if isinstance(arg, str) and arg.lstrip().startswith("{"):
        return json.loads(arg)
    with open(arg) as fh:
        return json.load(fh)


def run_panel(opts):
    arg = opts.get("panel") or opts.get("json")
    data = _load_json_arg(arg)
    if data is None:
        print("no input"); return 2

    values = _panel_from_json(data)
    age = opts.get("age")
    if age is None and data.get("age") is not None:
        age = float(data["age"])
    if age is None:
        at = None
        for key in ("measured_at", "date"):
            if data.get(key):
                at = _parse_date(str(data[key]))
        if at is None:
            at = date.today()
        if opts.get("birth_date"):
            age = age_at_date(opts["birth_date"], at)
        elif opts.get("birth_year") is not None:
            age = chronological_age(opts["birth_year"], at)
    if age is None:
        print("no age given (pass --age, --birth-year or an 'age' field)"); return 2

    missing = CRIT - set(values)
    res = {
        "complete": not missing,
        "chronological_age": round(age, 2),
        "phenoage": round(phenoage(values, age), 2) if not missing else None,
        "phenoage_acceleration": None,
        "markers_used": {k: values[k] for k in sorted(values)},
        "missing_markers": sorted(missing),
    }
    if not missing:
        res["phenoage_acceleration"] = round(res["phenoage"] - age, 2)
    else:
        res["note"] = "Cannot compute a valid PhenoAge until all 9 markers are present."

    if opts.get("format") == "json":
        print(json.dumps(res, indent=2))
    else:
        print_panel(res)
    return 0


def run_history(opts):
    path = opts["history"]
    birth_year = opts.get("birth_year")
    birth_date = opts.get("birth_date")

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        columns = {}
        for name in reader.fieldnames or []:
            kind, unit = parse_header(name)
            if kind:
                columns[name] = (kind, unit)
        rows = list(reader)

    if not any(k == "__date__" for k, _ in columns.values()):
        print("history CSV needs a date column (date or measured_at)"); return 2

    draws = []
    for row in rows:
        try:
            draws.append(_panel_from_row(row, columns, birth_year, birth_date))
        except ValueError as exc:
            print(f"row {row.get(next(iter(columns))) or '?'}: {exc}"); return 2

    draws.sort(key=lambda d: (d["measured_at"], d["id"]))
    changes = []
    for cur in draws:
        prev = previous_draw(draws, cur)
        if prev is not None:
            changes.append(explain_change(prev, cur))

    if opts.get("format") == "json":
        out = {
            "draws": [{k: v for k, v in d.items() if k != "date"} for d in draws],
            "changes": changes,
        }
        print(json.dumps(out, indent=2))
    else:
        print_trend(draws, changes)
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    opts = {"format": "text"}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a.startswith("--age="):
            opts["age"] = float(a.split("=", 1)[1]); i += 1
        elif a == "--age":
            opts["age"] = float(argv[i + 1]); i += 2
        elif a.startswith("--birth-year="):
            opts["birth_year"] = int(a.split("=", 1)[1]); i += 1
        elif a == "--birth-year":
            opts["birth_year"] = int(argv[i + 1]); i += 2
        elif a.startswith("--birth-date="):
            opts["birth_date"] = date.fromisoformat(a.split("=", 1)[1]); i += 1
        elif a == "--birth-date":
            opts["birth_date"] = date.fromisoformat(argv[i + 1]); i += 2
        elif a == "--history":
            opts["history"] = argv[i + 1]; i += 2
        elif a == "--format":
            opts["format"] = argv[i + 1]; i += 2
        elif a == "--json":
            nxt = argv[i + 1] if i + 1 < len(argv) else None
            if nxt and not nxt.startswith("--"):
                opts["json"] = nxt; i += 2
            else:
                opts["json"] = True; i += 1
        elif a.lstrip().startswith("{"):
            opts["panel"] = a; i += 1
        else:
            opts["panel"] = a; i += 1

    try:
        if opts.get("history"):
            return run_history(opts)
        return run_panel(opts)
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"error: {exc}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
