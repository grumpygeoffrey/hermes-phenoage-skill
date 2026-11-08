---
name: phenoage-calculator
description: Compute blood PhenoAge (Levine 2018 biological-age clock) and PhenoAgeAcceleration entirely locally from 9 standard clinical blood markers. No API key, no quota, no data leaves the machine. Use when asked for biological age, PhenoAge, an aging clock, or an age estimate from blood work.
version: 1.0.0
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

`age` is required too. Or pass `--age=52`.

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
`{ "age": 52, "markers": [ ... ] }`.

If markers are missing the script prints `"complete": false` with
`missing_markers` — never a made-up number.

## Caveats

- **Censored CRP** (lab prints `<0.3`) breaks the `ln`; a real CRP > 0 is required.
- Recompute only on **fasted, same-day** draws; don't mix draws across dates.
- This is for personal longevity research and education — not a diagnosis, not a
  substitute for a clinician. Outputs may be wrong. Not reviewed by any clinician.

## License & Attribution

MIT License. The verified implementation is derived from
[Insilico-org/plsi](https://github.com/Insilico-org/plsi) (`src/plsi/map/derive.py`),
MIT-licensed. The formula is public-domain science after Levine 2018
(doi:10.18632/aging.101414) and was cross-checked against the BioAge R reference
(Kwon & Belsky 2021, doi:10.1007/s11357-021-00480-5). Full MIT text in `LICENSE`.
