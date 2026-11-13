# hermes-phenoage-skill

A **Hermes skill** that computes the **Levine 2018 blood PhenoAge** — a validated
biological-age clock — from 9 standard clinical blood markers. Entirely local,
zero external dependencies, **no API key, no quota, and your blood data never
leaves your machine**.

v2.0 adds what the calculators usually lack: **trends across draws and a driver
decomposition** that says which marker moved your biological age, and by how many
years.

Contrast with the incumbent hub entry for the same function: a thin wrapper over a
paid third-party API. This skill runs the exact cross-verified formula locally.

## Layout

```
skills/phenoage-calculator/
├── SKILL.md            skill instructions & trigger routing
├── LICENSE              MIT (with upstream attribution)
└── scripts/
    └── phenoage.py       zero-dependency calculator, v2.0
tests/
└── test_phenoage.py     23 stdlib tests (maths, units, drivers, CLI)
```

## Install (for Hermes users)

```bash
hermes skills tap add grumpygeoffrey/hermes-phenoage-skill
hermes skills install phenoage-calculator
```

## Quick run

```bash
python3 skills/phenoage-calculator/scripts/phenoage.py \
  --birth-year 1974 markers.json
```

Prints `phenoage` and `phenoage_acceleration`. If any of the 9 markers is missing it
returns `"complete": false` with `missing_markers` — it will not fake a partial
number.

## Trend across draws

```bash
python3 skills/phenoage-calculator/scripts/phenoage.py --history draws.csv --birth-year 1974
```

```csv
date,birth_year,albumin_g_dl,creatinine_mg_dl,glucose_mg_dl,crp_mg_l,lymphocyte_pct,mcv_fl,rdw_pct,alp_u_l,wbc_10_3_ul
2025-11-20,1974,4.4,1.00,92,1.6,28.5,89.0,13.8,72,6.2
2026-11-13,1974,4.5,0.96,90,1.1,30.1,88.5,13.2,68,5.9
```

```
Blood PhenoAge — 2 draws

DATE            AGE  PHENOAGE    ACCEL
2025-11-20     51.4     48.11    -3.27
2026-11-13     52.4     44.77    -7.59

Change 2025-11-20 -> 2026-11-13: -3.34 y
  rdw            13.80 -> 13.20 %        -2.20 y
  age            51.4 -> 52.4            +0.87 y
  crp            0.16 -> 0.11 mg/dL      -0.40 y
```

Single-marker effects do not add exactly (the clock is non-linear in CRP); a
`combined` residual absorbs the remainder so the effects always sum to the total.

## Tests

```bash
python3 -m unittest discover -s tests -v
```

## Privacy / safety

- Runs byte-for-byte on the local machine. No network, no telemetry, no analytics.
- For personal longevity research and education; not a diagnosis and not a
  substitute for a clinician. Automated output is not physician-reviewed.

## License

MIT. Formula: Levine 2018 (doi:10.18632/aging.101414). Verified implementation
derived from Insilico-org/plsi (MIT); driver decomposition and UCUM unit scaling
adapted from brainarts-lab/hebe-longevity (MIT, Copyright (c) 2026 BrainArts).
See `LICENSE`.
