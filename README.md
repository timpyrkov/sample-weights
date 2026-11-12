# sample-weights: NHANES sample weights, expert review

One correlation in 9 cases. Each case is computed with the survey design (sample weights, strata, PSUs)
next to its naive twin (no weights). **For review by an expert in statistics (survey methods) and
bioinformatics: are the weights and the design handled correctly?**

## The variable pair

The two variables always have **different sample weights**:

| Variable | Measured in | Weight column (pynhanes) | NHANES name |
|---|---|---|---|
| BMI (kg/m2) | exam, everyone examined | `Sample weight (exam)` | WTMEC2YR |
| Glucose, fasting (mmol/L) | morning fasting subsample (~1/3) | `Sample weight (fasting)` | WTSAF2YR |

Analysis weight: **fasting** (the smaller subsample). Population: adults 20+, not pregnant.

## Rules

- Work in your own branch: `git checkout -b review/<your-name>`. Never commit to `master`.
- `data/` and `output/` are excluded from git sync (`.gitignore`). Do not add them.
- Code changes are welcome as proposals in your branch, each with its reason.

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
| 1 | One survey | 2017-18 | 0.196 / 0.0364 / 15 / 2,242 | 0.183 / 2,374 |
| 2 | Strata and PSUs ignored, weights kept: SE 0.0313, df 2,241 | 2017-18 | 0.196 / 0.0364 / 15 / 2,242 | 0.183 / 2,374 |
| 3 | Wrong weight: exam instead of fasting (r 0.202, n 2,374) | 2017-18 | 0.196 / 0.0364 / 15 / 2,242 | 0.183 / 2,374 |
| 4 | Two surveys, weight / 2 (undivided: same r, SE) | 2015-18 | 0.199 / 0.0219 / 30 / 4,471 | 0.175 / 4,773 |
| 5 | All 2-year surveys, weight / 8 | 2003-18 | 0.222 / 0.0102 / 124 / 17,620 | 0.183 / 18,991 |
| 6 | Pre vs post COVID, compared: 2017-18 | 2017-18 | 0.196 / 0.0364 / 15 / 2,242 | 0.183 / 2,374 |
|   | 2021-23 | 2021-23 | 0.163 / 0.0260 / 15 / 2,907 | 0.152 / 3,147 |
|   | difference | | 0.034 / 0.0447 / 30 | 0.031 |
| 7 | Pre + post COVID, pooled, weight / 2 | 2017-18, 2021-23 | 0.178 / 0.0222 / 30 / 5,149 | 0.166 / 5,521 |
| 8 | Filter age 30+ (rows deleted first: same) | 2017-18 | 0.199 / 0.0393 / 15 / 1,921 | 0.184 / 2,032 |
| 9 | Filter Mexican American (rows deleted first: SE 0.0264, df 11) | 2017-18 | 0.151 / 0.0280 / 15 / 325 | 0.142 / 334 |

## Review questions

1. **Weight choice.** Is the fasting weight right for BMI x fasting glucose (cases 1, 3)? Per NCHS guidelines.
2. **Strata and PSUs.** Is the Taylor-linearised SE with df = PSUs - strata right (case 2)?
3. **Rescaling.** Is dividing pooled weights by the number of surveys right, and does it leave r and SE unchanged (cases 4, 5, 7)?
4. **Survey choice.** Starting at 2003 (1999-2002 need 4-year weights); pooling 2021-23 with 2017-18 (cases 5-7).
5. **Independent check.** Do cases 1, 4, 9 agree with R `survey` / Stata `svy` / SAS (r, SE, df)?
6. **Filters.** Is a subgroup correctly analysed as a domain, not by deleting rows (cases 8, 9)?
7. **Inference.** Fisher-z CI with t(df) quantile, p from t(df): acceptable for a design-based r?
8. **Missing.** Anything a careful NHANES analysis of this pair should do and does not?

## Feedback format

Everything in your branch, as one pull request (do not merge). Suggested time: 2-4 hours, questions 1-5 first.

1. **`review/<your-name>.md`** (required, 1-2 pages, Markdown):
   - verdict table: question | cases | correct / incorrect / debatable | one-line reason | reference;
   - per incorrect or debatable point: what is wrong, effect on r / SE / p, suggested fix,
     severity (critical: changes conclusions; minor: changes numbers; style).
2. **Independent check** (optional, most valuable): e.g. an R `survey` script recomputing cases 1, 4, 9.
3. **Code changes** (optional, clear errors only): separate commits, reason in each message.
