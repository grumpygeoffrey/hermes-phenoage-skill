# hermes-phenoage-skill

A **Hermes skill** that computes the **Levine 2018 blood PhenoAge** — a validated
biological-age clock — from 9 standard clinical blood markers. Entirely local,
zero external dependencies, **no API key, no quota, and your blood data never
leaves your machine**.

Contrast with the incumbent hub entry for the same function: a thin wrapper over a
paid third-party API. This skill runs the exact cross-verified formula locally.

## Layout

```
skills/phenoage-calculator/
├── SKILL.md            skill instructions & trigger routing
├── LICENSE              MIT (with upstream attribution)
└── scripts/
    └── phenoage.py       zero-dependency calculator
```

## Install (for Hermes users)

```bash
hermes skills tap add grumpygeoffrey/hermes-phenoage-skill
hermes skills install phenoage-calculator
```

## Quick run

```bash
python3 skills/phenoage-calculator/scripts/phenoage.py \
  --age 52 --json='{"markers":[...]}'
```

Prints `phenoage` and `phenoage_acceleration`. If any of the 9 markers is missing it
returns `"complete": false` with `missing_markers` — it will not fake a partial
number.

See `skills/phenoage-calculator/SKILL.md` for the full input table and caveats.

## Privacy / safety

- Runs byte-for-byte on the local machine. No network, no telemetry, no analytics.
- For personal longevity research and education; not a diagnosis and not a
  substitute for a clinician. Automated output is not physician-reviewed.

## License

MIT. Formula: Levine 2018 (doi:10.18632/aging.101414). Verified implementation
derived from Insilico-org/plsi (MIT). See `LICENSE`.
