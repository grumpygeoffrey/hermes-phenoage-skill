---
name: phenoage-calculator
description: Compute blood PhenoAge (Levine 2018 biological-age clock) and PhenoAgeAcceleration entirely locally from 9 standard clinical blood markers. No API key, no quota, no data leaves the machine. Use when asked for biological age, PhenoAge, an aging clock, or an age estimate from blood work.
version: 2.0.0
author: grumpygeoffrey
license: MIT
platforms: [linux, macos]
tags: [longevity, biological-age, phenoage, aging-clock, blood-work, health]
related_skills: [blood-work-analyzer, longevity-health-analysis]
---

# PhenoAge Calculator

Compute the **Levine 2018 blood PhenoAge** — a validated biological-age clock — from
9 standard clinical blood markers plus chronological age. Runs as a zero-dependency
Python script, **entirely on the local machine**: no API, no key, no network, nothing
egresses.

This is the *blood* PhenoAge, not biolearn's DNAm (methylation) PhenoAge — no
epigenetic array required.

## When to use

- User asks "what's my biological age?" / "my PhenoAge" / "run an aging clock".
- User provides (or has in their data) the 9 blood markers below + age.
- Computing PhenoAgeAcceleration (= PhenoAge − chronological age) after a change in
  diet, retatrutide/tirzepatide, or lifestyle.

## When NOT to use

- Only some markers available → the script **refuses** to return a partial number
  (it lists the missing ones). Do not guess or fill from older draws.
- The user expects an overall biological-age composite spanning cognition/strength —
  PhenoAge measures mortality-related biology only.
- User expects GrimAge/DunedinPACE — those need a methylation array.

## Required inputs (canonical units)

| key | units | notes |
|---|---|---|
| `albumin` | g/L | |
| `creatinine` | µmol/L | or `creatinine_mgd` in mg/dL |
| `glucose` | mmol/L | or `glucose_mg` in mg/dL |
| `crp` | mg/dL | or `crp_mgl` in mg/L; must be > 0 (formula takes ln CRP) |
| `lymphocyte_pct` | % of WBC | |
| `mcv` | fL | |
| `rdw` | % (red cell distribution width) | |
| `alp` | U/L | |
| `wbc` | 10⁹/L | |

`age` is required too. Or pass `--age 52`. If you know the birth year/date and the
draw date, prefer `--birth-year 1974` / `--birth-date 1974-06-15` — the script then
computes the exact age at the draw (a 0.5-year error in age moves PhenoAge by ~0.04 y).

## Usage

```bash
python3 scripts/phenoage.py --age 52 '{"markers":[
  {"key":"albumin","value":45},{"key":"creatinine","value":84},{"key":"glucose","value":5.0},
  {"key":"crp","value":0.11},{"key":"lymphocyte_pct","value":30},{"key":"mcv","value":91},
  {"key":"rdw","value":12.5},{"key":"alp","value":70},{"key":"wbc","value":6.2}]}'
```

Output:
```json
{
  "complete": true,
  "chronological_age": 52.0,
  "phenoage": 42.78,
  "phenoage_acceleration": -9.22,
  "markers_used": { ... }
}
```

To run from a file: `python3 scripts/phenoage.py markers.json` where the JSON has
`{ "age": 52, "markers": [ ... ] }`. A flat `{"markers": {"crp": 0.11, ...}}` and an
optional per-marker `"unit"` field are also accepted (so labs can be copied verbatim,
e.g. `{"key":"glucose","value":92,"unit":"mg/dL"}`).

If markers are missing the script prints `"complete": false` with
`missing_markers` — never a made-up number.

### Trend across draws + drivers

```bash
python3 scripts/phenoage.py --history draws.csv --birth-year 1974
python3 scripts/phenoage.py --history draws.csv --format json
```

CSV header needs a date column (`date` or `measured_at`) and the 9 markers. Header
units are honoured: `albumin_g_dl`, `creatinine_mg_dl`, `glucose_mg_dl`, `crp_mg_l`,
`wbc_10_3_ul`, `alp_u_l`, `mcv_fl`, `rdw_pct`, `lymphocyte_pct` — or canonical names
for canonical units. A `birth_year`, `birth_date` or `age` column also works.

Output: a table of PhenoAge and acceleration per draw, then for each consecutive pair
a **driver decomposition** — which marker moved PhenoAge and by how many years:

```
Change 2025-11-20 -> 2026-11-13: -3.34 y
  rdw            13.80 -> 13.20 %        -2.20 y
  age            51.4 -> 52.4            +0.87 y
  crp            0.16 -> 0.11 mg/dL      -0.40 y
  ...
  combined       (non-linearity)          +0.00 y
```

Single-marker effects do not add up exactly (the clock is non-linear in CRP, and age
enters separately); the `combined` residual absorbs the remainder so the effects always
sum to the real total. `--format json` returns `draws` + `changes` for downstream use.

## Caveats

- **Censored values** (lab prints `<0.3`, `>5`) are rejected with an explicit error.
  Substituting the bound would silently bias the result, so get the number or omit it.
- Recompute only on **fasted, same-day** draws; don't mix draws across dates.
- This is for personal longevity research and education — not a diagnosis, not a
  substitute for a clinician. Outputs may be wrong. Not reviewed by any clinician.

## License & Attribution

MIT License. The verified implementation is derived from
[Insilico-org/plsi](https://github.com/Insilico-org/plsi) (`src/plsi/map/derive.py`),
MIT-licensed. Driver decomposition and UCUM unit scaling are adapted from
[brainarts-lab/hebe-longevity](https://github.com/brainarts-lab/hebe-longevity)
(MIT, Copyright (c) 2026 BrainArts), which vendors PLSI under the same terms. The formula is public-domain science after Levine 2018
(doi:10.18632/aging.101414) and was cross-checked against the BioAge R reference
(Kwon & Belsky 2021, doi:10.1007/s11357-021-00480-5). Full MIT text in `LICENSE`.
