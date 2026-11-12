# sample-weights: NHANES sample weights, QA task

One correlation in 8 cases. Each case is computed the correct way (sample weights, strata, PSUs)
next to its naive twin (no weights). **Task: tell us whether the weights are handled correctly.**

## The variable pair

The two variables always have **different sample weights**:

| Variable | Measured in | Weight column (pynhanes) | NHANES name |
|---|---|---|---|
| BMI (kg/m2) | exam, everyone examined | `Sample weight (exam)` | WTMEC2YR |
| Glucose, fasting (mmol/L) | morning fasting subsample (~1/3) | `Sample weight (fasting)` | WTSAF2YR |

Analysis weight: **fasting** (the smaller subsample). Population: adults 20+, not pregnant.

## Rules

- Work in your own branch: `git checkout -b qa/<your-name>`. Never commit to `master`.
- `data/` and `output/` are excluded from git sync (`.gitignore`). Do not add them.
- Do not edit the code to make a check pass. Report it.

## Run

Needs Python 3.12 and [uv](https://docs.astral.sh/uv/).

```bash
uv sync
uv run python fetch_data.py
uv run python correlations.py
```

- `fetch_data.py`: downloads one NHANES table (pynhanes release, ~30 MB) into `data/`, checks its sha256.
- `correlations.py`: runs all cases in seconds; `correlations.py 3 5` runs some.
- Output: results on screen; `output/case_N_*.png` (one figure per case); `output/results.csv`.
- Method details: docstrings in `correlations.py`.

## Cases (expected results)

| # | Case | Surveys | Weighted r / SE / df / n | Naive r / n |
|---|---|---|---|---|
| 1 | One survey (+ weights without strata/PSU: SE 0.0313) | 2017-18 | 0.196 / 0.0364 / 15 / 2,242 | 0.183 / 2,374 |
| 2 | Wrong weight: exam instead of fasting (r 0.202, n 2,374) | 2017-18 | 0.196 / 0.0364 / 15 / 2,242 | 0.183 / 2,374 |
| 3 | Two surveys, weight / 2 (undivided: same r, SE) | 2015-18 | 0.199 / 0.0219 / 30 / 4,471 | 0.175 / 4,773 |
| 4 | All 2-year surveys, weight / 8 | 2003-18 | 0.222 / 0.0102 / 124 / 17,620 | 0.183 / 18,991 |
| 5 | Pre vs post COVID, compared: 2017-18 | 2017-18 | 0.196 / 0.0364 / 15 / 2,242 | 0.183 / 2,374 |
|   | 2021-23 | 2021-23 | 0.163 / 0.0260 / 15 / 2,907 | 0.152 / 3,147 |
|   | difference | | 0.034 / 0.0447 / 30 | 0.031 |
| 6 | Pre + post COVID, pooled, weight / 2 | 2017-18, 2021-23 | 0.178 / 0.0222 / 30 / 5,149 | 0.166 / 5,521 |
| 7 | Filter age 60+ (rows deleted first: same) | 2017-18 | 0.172 / 0.0418 / 15 / 849 | 0.175 / 909 |
| 8 | Filter Mexican American (rows deleted first: SE 0.0264, df 11) | 2017-18 | 0.151 / 0.0280 / 15 / 325 | 0.142 / 334 |

## What to check

1. **Weight choice.** Fasting weight is right for BMI x fasting glucose (cases 1, 2). Check against NHANES docs.
2. **Rescaling.** Dividing pooled weights by the number of surveys is right, and does not change r or SE (cases 3, 4, 6).
3. **Survey choice.** Starting at 2003 (1999-2002 need 4-year weights); pooling 2021-23 with 2017-18 (cases 4-6).
4. **SE and df.** Recompute cases 1, 3, 8 with R `survey` / Stata `svy` / SAS. r, SE, df must agree.
5. **Filters.** Subgroup as a domain, not deleted rows (cases 7, 8).
6. **Naive twin.** Equals `scipy.stats.pearsonr` on the same rows.

## Reporting

Per finding: case, what you expected (with source), what you got. Commit notes to your branch, push, open a pull request. Do not merge.
