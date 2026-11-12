"""Survey-weighted correlation of BMI and fasting glucose in NHANES, next to its naive twin.

    uv run python correlations.py          # all 8 cases
    uv run python correlations.py 3 5      # some cases

The pair. The two variables come from different subsamples, so they have different weights:

    BMI (kg/m2)                 exam, everyone examined          Sample weight (exam)     WTMEC2YR
    Glucose, fasting (mmol/L)   morning fasting subsample, ~1/3  Sample weight (fasting)  WTSAF2YR

Rules every weighted analysis follows (NCHS analytic guidelines):

1. Weight. Use the weight of the smallest subsample among the variables: the fasting weight.
   Only it corrects for who was picked for the morning session and came fasting.
2. Pooling surveys. Each survey keeps its own 2-year weight, divided by the number of surveys.
   The division (any constant rescaling) does not change r or its SE; it only keeps the weights
   summing to the population. 1999-2000 + 2001-2002 would need 4-year weights: not used here.
3. Standard error. NHANES samples areas (PSUs) within strata, then people. People of one PSU
   resemble each other, so the SE comes from PSU totals within strata (Taylor linearisation),
   with df = PSUs - strata (15 per survey), not n - 2.
4. Filters are domains. Rows outside the subgroup get zero weight but stay in the design, so no
   PSU or stratum is lost. Deleting rows before computing the SE is a mistake.

The naive twin ignores all four: no weights, rows independent, df = n - 2 (scipy.stats.pearsonr).
It also keeps rows with a glucose value but zero fasting weight, so its n is larger.

Population in every case: adults aged 20+, pregnant women left out.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # no display: figures go straight to files
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import fetch_data

OUTPUT_DIR = Path(__file__).parent / "output"
X, Y = "BMI (kg/m2)", "Glucose, fasting (mmol/L)"
EXAM, FASTING = "Sample weight (exam)", "Sample weight (fasting)"
SURVEY, STRATA, PSU = "Survey", "Survey strata", "Survey PSU"
AGE, PREGNANT, ETHNICITY = "Age", "Pregnancy status", "Ethnicity"  # Ethnicity 1: Mexican American
COLUMNS = [SURVEY, STRATA, PSU, AGE, PREGNANT, ETHNICITY, X, Y, EXAM, FASTING]


# ---- data and weights ------------------------------------------------------------------

@lru_cache(maxsize=1)
def _table() -> pd.DataFrame:
    path = fetch_data.DATA_DIR / fetch_data.TABLE
    with gzip.open(path, "rt") as f:  # line 1 categories, line 2 names, line 3 "SEQN"
        f.readline()
        names = f.readline().rstrip("\n").split(";")
    names[0] = "SEQN"
    return pd.read_csv(path, sep=";", skiprows=3, header=None, names=names, usecols=COLUMNS)


def load(surveys: list[str]) -> pd.DataFrame:
    """Every participant of these surveys, all ages: the full design, not only the analysed rows."""
    table = _table()
    return table[table[SURVEY].isin(surveys)].reset_index(drop=True)


def weight(data: pd.DataFrame, column: str = FASTING, divide: bool = True) -> np.ndarray:
    """One weight per row: the survey's own column (0 if missing), / number of pooled surveys (rule 2)."""
    return data[column].fillna(0).to_numpy() / (data[SURVEY].nunique() if divide else 1)


def adults(data: pd.DataFrame) -> np.ndarray:
    """The analysed population: aged 20+ and not pregnant."""
    return ((data[AGE] >= 20) & (data[PREGNANT] != 1)).to_numpy()


# ---- estimators ------------------------------------------------------------------------

@dataclass
class Result:
    label: str      # "weighted", "naive", or a variant such as "mistake"
    how: str        # weight and variance used
    r: float
    se: float
    df: int
    n: int
    p: float
    lo: float       # 95% CI
    hi: float
    group: str = "All"
    xy: tuple | None = None  # (x, y, w) of the analysed rows, for the plot; None for a difference


def design_se(infl: np.ndarray, strata: np.ndarray, psu: np.ndarray) -> tuple[float, int]:
    """SE of an estimate from its per-row influence, and df = PSUs - strata (rule 3).

    Influence is summed per PSU. In a stratum with n_h PSUs the variance share is
    n_h / (n_h - 1) * sum (PSU total - stratum mean)^2. A stratum with one PSU adds nothing.
    """
    totals = pd.DataFrame({"h": strata, "p": psu, "z": infl}).groupby(["h", "p"])["z"].sum()
    by_h = totals.groupby(level="h")
    n_h = by_h.transform("size")
    keep = n_h > 1
    dev = (totals - by_h.transform("mean"))[keep]
    var = float((n_h[keep] / (n_h[keep] - 1) * dev ** 2).sum())
    return float(np.sqrt(var)), int(keep.sum() - (by_h.size() > 1).sum())


def weighted_corr(data: pd.DataFrame, w: np.ndarray, domain: np.ndarray, how: str,
                  label: str = "weighted", group: str = "All", design: bool = True) -> Result:
    """Weighted Pearson r over the domain, with the design-based SE (rules 1-4).

    `data` holds every row of the surveys; `domain` marks the analysed ones. Rows outside it,
    without both values or with zero weight get zero influence but keep their PSU in the SE.
    design=False: weights kept, but every row is its own PSU in one stratum (rule 3 ignored).
    """
    ok = domain & (w > 0) & data[[X, Y]].notna().all(axis=1).to_numpy()
    wk = np.where(ok, w, 0.0)
    x, y = data[X].fillna(0).to_numpy(), data[Y].fillna(0).to_numpy()
    W = wk.sum()
    mx, my = (wk * x).sum() / W, (wk * y).sum() / W
    zx = np.where(ok, (x - mx) / np.sqrt((wk * (x - mx) ** 2).sum() / W), 0)
    zy = np.where(ok, (y - my) / np.sqrt((wk * (y - my) ** 2).sum() / W), 0)
    r = float((wk * zx * zy).sum() / W)
    infl = wk * (zx * zy - 0.5 * r * (zx ** 2 + zy ** 2)) / W  # how much each row moves r
    n = int(ok.sum())
    if design:
        se, df = design_se(infl, data[STRATA].to_numpy(), data[PSU].to_numpy())
        how += ", strata + PSU"
    else:
        se, df = design_se(infl, np.zeros(len(data)), np.arange(len(data)))[0], n - 1
        how += ", rows independent"
    q = stats.t.ppf(0.975, df)
    z, se_z = np.arctanh(r), se / (1 - r ** 2)  # Fisher-z interval with the t quantile
    return Result(label, how, r, se, df, n, float(2 * stats.t.sf(abs(r / se), df)),
                  float(np.tanh(z - q * se_z)), float(np.tanh(z + q * se_z)), group,
                  (x[ok], y[ok], wk[ok]))


def naive_corr(data: pd.DataFrame, domain: np.ndarray, group: str = "All") -> Result:
    """The naive twin: scipy.stats.pearsonr on every row with both values. No weights, no design."""
    ok = domain & data[[X, Y]].notna().all(axis=1).to_numpy()
    x, y = data.loc[ok, X].to_numpy(), data.loc[ok, Y].to_numpy()
    res = stats.pearsonr(x, y)
    ci = res.confidence_interval()
    n = len(x)
    se = float(np.sqrt((1 - res.statistic ** 2) / (n - 2)))
    return Result("naive", "no weights, rows independent", float(res.statistic), se, n - 2, n,
                  float(res.pvalue), float(ci.low), float(ci.high), group, (x, y, np.ones(n)))


def difference(a: Result, b: Result) -> Result:
    """r of a minus r of b, for two independent samples (two surveys): variances and df add."""
    d, se, df = a.r - b.r, float(np.hypot(a.se, b.se)), a.df + b.df
    q = stats.t.ppf(0.975, df)
    return Result(a.label, a.how, d, se, df, a.n + b.n, float(2 * stats.t.sf(abs(d / se), df)),
                  d - q * se, d + q * se, f"{a.group} minus {b.group}")


# ---- output ----------------------------------------------------------------------------

RESULTS: list[dict] = []
STYLE = {"weighted": ("C0", "-"), "naive": ("C3", "-"), "mistake": ("C1", "--"),
         "weights only": ("C2", ":"), "not divided": ("C4", "--")}  # colour, dash


def _line(res: Result) -> str:
    name = "difference in r" if res.xy is None else "r"
    return (f"{name} = {res.r:+.3f} (95% CI {res.lo:+.3f} to {res.hi:+.3f}), SE = {res.se:.4f}, "
            f"p = {res.p:.2g}, n = {res.n:,}, df = {res.df:,}")


def show(case: int, name: str, title: str, surveys: list[str], panels: list[list[Result]],
         footer: tuple[Result, ...] = (), note: str = "") -> None:
    """Print a case, add it to results.csv and save output/case_N_name.png (one panel per list)."""
    print(f"\nCase {case}: {title}  [{', '.join(surveys)}]")
    for res in [r for panel in panels for r in panel] + list(footer):
        print(f"  {res.label:<12} {res.group}: {res.how}\n      {_line(res)}")
        RESULTS.append({"case": case, "name": name, "surveys": " ".join(surveys), "group": res.group,
                        "analysis": res.label, "how": res.how, "n": res.n, "r": res.r, "se": res.se,
                        "df": res.df, "ci_low": res.lo, "ci_high": res.hi, "p": res.p})
    if note:
        print(f"  Note: {note}")

    fig, axes = plt.subplots(1, len(panels), figsize=(7.5 * len(panels), 7.6), squeeze=False)
    rng = np.random.default_rng(0)
    for ax, panel in zip(axes[0], panels):
        x, y, _ = panel[0].xy
        rows = rng.choice(len(x), min(3000, len(x)), replace=False)
        ax.scatter(x[rows], y[rows], s=5, alpha=0.25, color="0.45", linewidths=0)
        grid = np.array([x.min(), x.max()])
        for res in panel:  # weighted least-squares line of each analysis
            px, py, pw = res.xy
            mx, my = np.average(px, weights=pw), np.average(py, weights=pw)
            slope = np.average((px - mx) * (py - my), weights=pw) / np.average((px - mx) ** 2, weights=pw)
            colour, dash = STYLE[res.label]
            ax.plot(grid, my + slope * (grid - mx), color=colour, ls=dash, lw=2, label=res.label)
        ax.set(xlabel=X, ylabel=Y, title=panel[0].group if len(panels) > 1 else "")
        ax.legend(loc="upper left", frameon=False)
        text = "\n".join(f"{r.label}: {r.how}\n  {_line(r)}" for r in panel)
        ax.text(0, -0.13, text, transform=ax.transAxes, va="top", family="monospace", fontsize=8.5)
    foot = [f"{r.label}, {r.group}: {_line(r)}" for r in footer] + ([f"Note: {note}"] if note else [])
    fig.text(0.02, 0.012, "\n".join(foot), family="monospace", fontsize=8.5, va="bottom")
    fig.suptitle(f"Case {case}: {title}\n{', '.join(surveys)}; points: up to 3,000 analysed rows",
                 fontsize=11)
    fig.subplots_adjust(left=0.09 if len(panels) == 1 else 0.06, right=0.98, top=0.9, bottom=0.3,
                        wspace=0.18)
    OUTPUT_DIR.mkdir(exist_ok=True)
    fig.savefig(OUTPUT_DIR / f"case_{case}_{name}.png", dpi=130)
    plt.close(fig)


# ---- the cases -------------------------------------------------------------------------

def case_1_one_survey() -> None:
    """Base case: one survey, the fasting weight, the design-based SE.

    A third line, "weights only", keeps the weights but treats rows as independent: r is the
    same (r depends on the weights only), the SE is too small and df is n - 1, not 15.
    Weights fix the estimate; strata and PSUs fix its uncertainty.
    """
    surveys = ["2017-2018"]
    data = load(surveys)
    w = weight(data)
    show(1, "one_survey", "one survey", surveys,
         [[weighted_corr(data, w, adults(data), "fasting weight"),
           weighted_corr(data, w, adults(data), "fasting weight", "weights only", design=False),
           naive_corr(data, adults(data))]])


def case_2_wrong_weight() -> None:
    """Case 1 with a mistake: the exam weight (the weight of BMI) instead of the fasting weight.

    The exam weight is positive for everyone examined, so it also keeps people outside the
    fasting subsample (n grows) and does not correct for who came fasting.
    """
    surveys = ["2017-2018"]
    data = load(surveys)
    show(2, "wrong_weight", "the wrong weight (exam instead of fasting)", surveys,
         [[weighted_corr(data, weight(data), adults(data), "fasting weight"),
           weighted_corr(data, weight(data, EXAM), adults(data), "exam weight", "mistake"),
           naive_corr(data, adults(data))]])


def case_3_two_surveys() -> None:
    """Two adjacent surveys pooled: each keeps its own fasting weight, divided by 2.

    "not divided" uses the raw weights: r and SE are identical, because a constant factor
    cancels (rule 2). The division only matters for totals. Strata numbers differ between
    surveys, so df doubles to 30.
    """
    surveys = ["2015-2016", "2017-2018"]
    data = load(surveys)
    show(3, "two_surveys", "two surveys pooled (weight / 2)", surveys,
         [[weighted_corr(data, weight(data), adults(data), "fasting weight / 2"),
           weighted_corr(data, weight(data, divide=False), adults(data), "fasting weight",
                         "not divided"),
           naive_corr(data, adults(data))]])


def case_4_all_surveys() -> None:
    """All eight 2-year surveys 2003-2018 pooled, weight / 8.

    Starts at 2003: 1999-2000 + 2001-2002 need the special 4-year weights. The naive SE keeps
    shrinking with n; the design SE is limited by the number of PSUs (df = 124, not n - 2).
    """
    surveys = [f"{y}-{y + 1}" for y in range(2003, 2018, 2)]
    data = load(surveys)
    show(4, "all_surveys", "all 2-year surveys 2003-2018 pooled (weight / 8)", surveys,
         [[weighted_corr(data, weight(data), adults(data), "fasting weight / 8"),
           naive_corr(data, adults(data))]])


def case_5_pre_post_covid() -> None:
    """Before (2017-2018) and after (2021-2023) the pandemic, compared, not pooled.

    Each survey on its own with its own undivided weight (df 15 each). The surveys are
    independent samples, so the variance of the difference is the sum of the two variances.
    """
    surveys = ["2017-2018", "2021-2023"]
    wtd, nai = [], []
    for s in surveys:
        data = load([s])
        wtd.append(weighted_corr(data, weight(data), adults(data), "fasting weight", group=s))
        nai.append(naive_corr(data, adults(data), group=s))
    show(5, "pre_post_covid", "before and after the pandemic, compared", surveys,
         [[wtd[0], nai[0]], [wtd[1], nai[1]]], footer=(difference(*wtd), difference(*nai)))


def case_6_across_covid() -> None:
    """2017-2018 and 2021-2023 pooled, weight / 2, as in case 3.

    The arithmetic is the same, but 2021-2023 was sampled and weighted differently after the
    pandemic and 2019-2020 is missing: the result averages two periods, it is not a 2017-2023
    estimate. (The 2017-March 2020 file is not in the data; it would replace 2017-2018.)
    """
    surveys = ["2017-2018", "2021-2023"]
    data = load(surveys)
    show(6, "across_covid", "pooled across the pandemic (weight / 2)", surveys,
         [[weighted_corr(data, weight(data), adults(data), "fasting weight / 2"),
           naive_corr(data, adults(data))]],
         note="two periods and two designs averaged; compare with case 5")


def _filtered(case: int, name: str, title: str, group: str, member) -> None:
    """Cases 7 and 8: a filter as a domain (correct) and with the other rows deleted first (mistake)."""
    surveys = ["2017-2018"]
    data = load(surveys)
    w = weight(data)
    domain = adults(data) & member(data).to_numpy()
    kept = data[domain].reset_index(drop=True)
    show(case, name, title, surveys,
         [[weighted_corr(data, w, domain, "fasting weight", group=group),
           weighted_corr(kept, w[domain], np.ones(len(kept), bool),
                         "fasting weight, other rows deleted", "mistake", group),
           naive_corr(data, domain, group)]])


def case_7_age_filter() -> None:
    """Filter by age (60+), analysed as a domain (rule 4).

    "mistake" deletes the other rows first. Older adults live in every PSU, so here no PSU is
    lost and both give the same SE and df: harmless by luck, see case 8.
    """
    _filtered(7, "age_filter", "filtered by age (60+)", "Age 60+", lambda d: d[AGE] >= 60)


def case_8_ethnicity_filter() -> None:
    """Filter by ethnicity (Mexican American), analysed as a domain (rule 4).

    Mexican Americans are absent from some PSUs. Deleting the other rows first drops those
    PSUs, strata left with one PSU add no variance: the SE shrinks and df falls from 15 to 11.
    """
    _filtered(8, "ethnicity_filter", "filtered by ethnicity (Mexican American)",
              "Mexican American", lambda d: d[ETHNICITY] == 1)


CASES = [case_1_one_survey, case_2_wrong_weight, case_3_two_surveys, case_4_all_surveys,
         case_5_pre_post_covid, case_6_across_covid, case_7_age_filter, case_8_ethnicity_filter]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="*", type=int, help="case numbers 1-8 (default: all)")
    args = ap.parse_args()
    bad = fetch_data.problems()
    if bad:
        sys.exit(f"data/ is not ready ({'; '.join(bad)}). Run: uv run python fetch_data.py")
    for i in args.cases or range(1, len(CASES) + 1):
        CASES[i - 1]()
    with open(OUTPUT_DIR / "results.csv", "w", newline="") as f:
        out = csv.DictWriter(f, fieldnames=list(RESULTS[0]))
        out.writeheader()
        out.writerows(RESULTS)
    print(f"\nFigures and results.csv are in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
