# sample-weights: weighted correlations in NHANES, QA task

One correlation, BMI x fasting glucose, in 16 situations. Each is computed the correct way
(sample weights, strata, PSUs) and the naive way (no weights). The task: check that the correct way is correct.

Files: `correlations.py` (all the analysis; the explanations are in its docstrings), `fetch_data.py` (gets the data).

## Rules

- Work in your own branch: `git checkout -b qa/<your-name>`. Do not commit to `master`.
- `data/` and `output/` are excluded from git (`.gitignore`). Do not add them.
- Do not edit `correlations.py` to make a check pass. Report the finding.

## Setup and run

Needs Python 3.12 and [uv](https://docs.astral.sh/uv/).

```bash
uv sync
uv run python fetch_data.py
uv run python correlations.py
```

- `uv sync` installs the dependencies, including pynhanes 1.1.0.
- `fetch_data.py` creates `data/` (2 public NHANES tables from the pynhanes release, about 30 MB) and verifies their sha256.
- `correlations.py` runs all cases (under a minute). `correlations.py 3 7` runs some; `--list` lists them.
- Output: results on screen, one figure per case in `output/case_NN_*.png`, all numbers in `output/results.csv`.

## Cases

Adults 20+, not pregnant. 2017-2018 unless stated. Weighted: r, SE, df, n. Naive: r, n.

| # | Case | Weighted r | SE | df | n | Naive r | n |
|---|---|---|---|---|---|---|---|
| 1 | One survey, fasting weight | 0.196 | 0.0364 | 15 | 2,242 | 0.183 | 2,374 |
| 2 | Mistake: exam weight (r 0.202, n 2,374) | 0.196 | 0.0364 | 15 | 2,242 | 0.183 | 2,374 |
| 3 | Three weight levels: age x sedentary time, interview weight | 0.050 | 0.0241 | 15 | 5,470 | 0.039 | 5,470 |
|   | age x BMI, exam weight | 0.025 | 0.0230 | 15 | 5,120 | -0.004 | 5,120 |
|   | BMI x glucose, fasting weight | 0.196 | 0.0364 | 15 | 2,242 | 0.183 | 2,374 |
| 4 | Weights without strata and PSUs (SE 0.0313, df 2,241) | 0.196 | 0.0364 | 15 | 2,242 | 0.183 | 2,374 |
| 5 | Weights rescaled to mean 1 (r, SE, df unchanged) | 0.196 | 0.0364 | 15 | 2,242 | 0.183 | 2,374 |
| 6 | Two surveys pooled, 2015-2018 | 0.199 | 0.0219 | 30 | 4,471 | 0.175 | 4,773 |
| 7 | Seven surveys pooled, 2005-2018 | 0.222 | 0.0110 | 109 | 15,766 | 0.184 | 16,977 |
| 8 | Before / after the pandemic: 2017-2018 | 0.196 | 0.0364 | 15 | 2,242 | 0.183 | 2,374 |
|   | 2021-2023 | 0.163 | 0.0260 | 15 | 2,907 | 0.152 | 3,147 |
|   | difference | 0.034 | 0.0447 | 30 | | 0.031 | |
| 9 | 2017-2018 + 2021-2023 pooled | 0.178 | 0.0222 | 30 | 5,149 | 0.166 | 5,521 |
| 10 | BMI x HbA1c, weight changes between surveys (mistake, exam weight: r 0.207) | 0.205 | 0.0123 | 30 | 10,542 | 0.185 | 10,542 |
| 11 | Subgroup: Mexican American (mistake, rows deleted first: SE 0.0264, df 11) | 0.151 | 0.0280 | 15 | 325 | 0.142 | 334 |
| 12 | Subgroup: age 40-59 (rows deleted first: same SE and df) | 0.207 | 0.0680 | 15 | 729 | 0.186 | 763 |
| 13 | By gender: men | 0.194 | 0.0510 | 15 | 1,077 | 0.215 | 1,151 |
|   | women | 0.223 | 0.0411 | 15 | 1,165 | 0.173 | 1,223 |
|   | difference | -0.030 | 0.0612 | 15 | | 0.043 | |
| 14 | Adjusted for age and gender | 0.207 | 0.0340 | 15 | 2,242 | 0.198 | 2,374 |
| 15 | Spearman | 0.317 | 0.0341 | 15 | 2,242 | 0.269 | 2,374 |
| 16 | Refused: 1999-2000 + 2001-2002 (needs 4-year weights); 2017-2020 (not in the data) | | | | | | |

## What to check

1. **Reproduction.** Your run gives the numbers above and the same `results.csv` on a second run.
2. **Weight choice.** Each case uses the weight NCHS prescribes for its variables and surveys (cases 1, 2, 3, 10).
   Compare with the NHANES documentation, not with the script.
3. **Weights vs design.** Case 4: r equals case 1, SE is smaller. Case 5: rescaling changes nothing.
4. **Pooling.** Each survey keeps its own weight, divided by the number of surveys (cases 6, 7, 9, 10).
   Is the claim "the division does not change r" true?
5. **Standard error and df.** Recompute at least cases 1, 6 and 11 independently
   (R `survey`, Stata `svy`, SAS, or your own code). r, SE, df, CI, p must agree.
6. **Subgroups.** Case 11: the domain result is correct, the rows-deleted result is wrong.
   Case 12: both agree; is the stated reason right? Case 13: the SE of the difference.
7. **Naive twin.** Equals `scipy.stats.pearsonr` / `spearmanr` on the same rows.
8. **Adjusted and Spearman** (cases 14, 15): method and SE are acceptable; state the limits you find.
9. **Across the pandemic** (cases 8, 9): the treatment of 2021-2023 and its caveat are right.
10. **Refusals** (case 16): both fire. Find other combinations that should be refused and are not.
11. **Figures.** The numbers on each figure equal `results.csv`; lines and labels are right.
12. **Docstrings.** Each says what the code does. Report any statement that is wrong or unclear.

## Reporting

One entry per finding: case number, what you expected (with source), what you got, how to reproduce.
Commit your notes and any checking scripts to your branch, push the branch, open a pull request. Do not merge.
