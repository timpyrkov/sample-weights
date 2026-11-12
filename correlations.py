"""Survey-weighted correlations in NHANES, each next to its naive unweighted twin.

    uv run python correlations.py           # all cases
    uv run python correlations.py 3 7       # some cases
    uv run python correlations.py --list    # the case list

One pair of variables measured on different subsamples:

    BMI (kg/m2)                 body measures, exam weight (everyone examined)
    Glucose, fasting (mmol/L)   morning fasting subsample (about 1/3), fasting weight

Each `case_NN_*` function below is one situation a researcher meets. Each prints its results,
adds rows to output/results.csv and saves output/case_NN_*.png.

The rules every case follows (NCHS analytic guidelines):

1. Weight. Use the weight of the smallest subsample among the variables: here the fasting
   weight. It corrects for unequal selection, non-response and, for subsamples, for who was
   picked and came fasting. `pick_weight` reads the weight from the pynhanes weights dictionary.
2. Pooling surveys. Each survey keeps its own 2-year weight, divided by the number of surveys.
   The division does not change r (a constant factor cancels); it keeps the weights summing
   to the average population. 1999-2000 with 2001-2002 needs special 4-year weights.
3. Standard error. NHANES samples areas (PSUs) within strata, then people. People of one PSU
   resemble each other, so the SE comes from PSU totals within strata (Taylor linearisation),
   and the degrees of freedom are PSUs - strata (15 per survey), not n - 2.
4. Subgroups are domains. Rows outside the subgroup get zero influence but stay in the
   design, so no PSU or stratum is lost. Do not delete rows before computing the SE.

The naive twin ignores all four: no weights, every row independent, df = n - 2. It also
keeps the rows whose weight is zero (a glucose value, but outside the fasting subsample),
so its n is larger. Case 4 adds a half-way analysis (weights, but no strata and PSUs) to
show which part of the difference comes from the weights and which from the design.

Population in every case: adults aged 20+, pregnant women left out.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # no display needed: figures go straight to files
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import fetch_data

DATA_DIR = fetch_data.DATA_DIR
OUTPUT_DIR = Path(__file__).parent / "output"

X, Y = "BMI (kg/m2)", "Glucose, fasting (mmol/L)"
HBA1C = "HbA1c (mmol/mol)"
SURVEY, STRATA, PSU = "Survey", "Survey strata", "Survey PSU"
AGE, GENDER, PREGNANT = "Age", "Gender", "Pregnancy status"  # Gender: 1 male, 0 female
ETHNICITY = "Ethnicity"  # 1 Mexican American
SEDENTARY = "Physical activity sedentary (min/day)"  # questionnaire, used in case 3 only
FASTING, BLOOD, EXAM, INTERVIEW = ("Sample weight (fasting)", "Sample weight (blood draw)",
                                   "Sample weight (exam)", "Sample weight (interview)")
SMALLEST_FIRST = [FASTING, BLOOD, EXAM, INTERVIEW]  # subsample size, smallest first
COLUMNS = [SURVEY, STRATA, PSU, AGE, GENDER, PREGNANT, ETHNICITY, SEDENTARY, X, Y, HBA1C,
           *SMALLEST_FIRST]
CI = 0.95


class Refused(Exception):
    """The analysis cannot be done correctly with this data; the message says why."""


# ---- data, weights, population ---------------------------------------------------------

@lru_cache(maxsize=1)
def _table() -> pd.DataFrame:
    path = DATA_DIR / "nhanes_userdata.csv.gz"
    with gzip.open(path, "rt") as f:  # line 1 categories, line 2 names, line 3 "SEQN"
        f.readline()
        names = f.readline().rstrip("\n").split(";")
    names[0] = "SEQN"
    return pd.read_csv(path, sep=";", skiprows=3, header=None, names=names, usecols=COLUMNS)


@lru_cache(maxsize=1)
def _weights_dict() -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / "nhanes_weights_dict.csv.gz", sep=";", index_col=0)


def load(surveys: list[str]) -> pd.DataFrame:
    """Every participant of these surveys, all ages: the full design, not only the analysed rows."""
    table = _table()
    missing = [s for s in surveys if s not in set(table[SURVEY])]
    if missing:
        raise Refused(f"no rows for {', '.join(missing)} in this snapshot")
    return table[table[SURVEY].isin(surveys)].reset_index(drop=True)


def pick_weight(variables: list[str], survey: str) -> str:
    """The weight column for these variables in one survey: that of the smallest subsample."""
    book = _weights_dict()
    used = []
    for v in variables:
        w = book.at[v, survey]
        if not isinstance(w, str):
            raise Refused(f"{v} was not measured in {survey}")
        used.append(w)
    return min(used, key=SMALLEST_FIRST.index)


def analysis_weight(data: pd.DataFrame, variables: list[str], force: str | None = None,
                    divide: bool = True) -> tuple[np.ndarray, str]:
    """One weight per row, and its description.

    Per survey: the column from `pick_weight` (or `force`, to show a mistake), divided by the
    number of surveys in `data` when they are pooled into one estimate (`divide`).
    """
    surveys = sorted(data[SURVEY].unique())
    if {"1999-2000", "2001-2002"} <= set(surveys):
        raise Refused("1999-2000 and 2001-2002 can only be pooled with the 4-year weights "
                      "(WTMEC4YR, WTSAF4YR); this snapshot has only 2-year weights")
    w = np.zeros(len(data))
    cols = {}
    for s in surveys:
        cols[s] = force or pick_weight(variables, s)
        rows = (data[SURVEY] == s).to_numpy()
        w[rows] = data.loc[rows, cols[s]].fillna(0).to_numpy()
    if divide:
        w = w / len(surveys)
    short = {s: c[len("Sample weight ("):-1] for s, c in cols.items()}
    text = (next(iter(short.values())) if len(set(short.values())) == 1
            else ", ".join(f"{c} in {s}" for s, c in short.items()))
    text += " weight"
    if divide and len(surveys) > 1:
        text += f" / {len(surveys)}"
    return w, text


def adults(data: pd.DataFrame) -> np.ndarray:
    """The analysed population: aged 20+ and not pregnant."""
    return ((data[AGE] >= 20) & ~(data[PREGNANT] == 1)).to_numpy()


# ---- estimators ------------------------------------------------------------------------

@dataclass
class Result:
    label: str          # "weighted", "naive", "weights only", "rescaled" or "mistake"
    how: str            # weight and variance used
    r: float
    se: float
    df: int
    n: int
    group: str = "All"
    x: str = ""         # names of the two variables
    y: str = ""
    # kept for the plot and for differences between groups
    xs: np.ndarray = field(default=None, repr=False)
    ys: np.ndarray = field(default=None, repr=False)
    w: np.ndarray = field(default=None, repr=False)
    ok: np.ndarray = field(default=None, repr=False)
    infl: np.ndarray = field(default=None, repr=False)
    is_difference: bool = False
    textbook_bounds: tuple[float, float] | None = None  # set by the naive analysis

    @property
    def p(self) -> float:
        return float(2 * stats.t.sf(abs(self.r / self.se), self.df))

    @property
    def bounds(self) -> tuple[float, float]:
        """Fisher-z interval with the t quantile (plain interval for a difference)."""
        if self.textbook_bounds:
            return self.textbook_bounds
        q = stats.t.ppf(0.5 + CI / 2, self.df)
        if self.is_difference:
            return self.r - q * self.se, self.r + q * self.se
        z, se_z = np.arctanh(self.r), self.se / (1 - self.r ** 2)
        return float(np.tanh(z - q * se_z)), float(np.tanh(z + q * se_z))


def design_variance(infl: np.ndarray, strata: np.ndarray, psu: np.ndarray) -> tuple[float, int]:
    """Variance of an estimate from its influence rows, and the design degrees of freedom.

    Influence rows are summed per PSU. Within each stratum with n_h PSUs the share of the
    variance is n_h / (n_h - 1) * sum (PSU total - stratum mean)^2. df = PSUs - strata.
    """
    totals = pd.DataFrame({"s": strata, "p": psu, "z": infl}).groupby(["s", "p"])["z"].sum()
    by_stratum = totals.groupby(level="s")
    n_h = by_stratum.transform("size")
    centred = totals - by_stratum.transform("mean")
    many = n_h > 1  # a stratum with a single PSU cannot show variation
    var = float((n_h[many] / (n_h[many] - 1) * centred[many] ** 2).sum())
    return var, int(many.sum() - by_stratum.ngroups + (n_h.groupby(level="s").first() == 1).sum())


def ridit(v: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted rank in (0, 1): the population share below each value plus half of its ties."""
    _, inv = np.unique(v, return_inverse=True)
    tot = np.bincount(inv, weights=w)
    below = np.concatenate([[0.0], np.cumsum(tot)[:-1]])
    return (below[inv] + 0.5 * tot[inv]) / w.sum()


def residuals(v: np.ndarray, covariates: np.ndarray, w: np.ndarray) -> np.ndarray:
    """What is left of v after a weighted least-squares fit on the covariates."""
    C = np.column_stack([np.ones(len(v)), covariates]) * np.sqrt(w)[:, None]
    beta, *_ = np.linalg.lstsq(C, v * np.sqrt(w), rcond=None)
    return v - np.column_stack([np.ones(len(v)), covariates]) @ beta


def _corr(data, x, y, w, domain, method, adjust, strata, psu, label, how, group) -> Result:
    """Shared by `weighted_corr` and `naive_corr`: r and its linearised SE for given weights."""
    cols = [x, y, *adjust]
    ok = domain & (w > 0) & data[cols].notna().all(axis=1).to_numpy()
    xs, ys = np.zeros(len(data)), np.zeros(len(data))
    xs[ok], ys[ok] = data.loc[ok, x], data.loc[ok, y]
    wk = np.where(ok, w, 0.0)
    if method == "spearman":
        xs[ok], ys[ok] = ridit(xs[ok], wk[ok]), ridit(ys[ok], wk[ok])
    if adjust:
        cov = data.loc[ok, adjust].to_numpy(float)
        xs[ok], ys[ok] = residuals(xs[ok], cov, wk[ok]), residuals(ys[ok], cov, wk[ok])
    W = wk.sum()
    mx, my = (wk * xs).sum() / W, (wk * ys).sum() / W
    sx = np.sqrt((wk * (xs - mx) ** 2).sum() / W)
    sy = np.sqrt((wk * (ys - my) ** 2).sum() / W)
    zx, zy = (xs - mx) / sx, (ys - my) / sy
    r = float((wk * zx * zy).sum() / W)
    # how much each row pulls r (the ranks and residuals are treated as fixed)
    infl = wk * (zx * zy - 0.5 * r * (zx ** 2 + zy ** 2)) / W
    var, df = design_variance(infl, strata, psu)
    return Result(label, how, r, float(np.sqrt(var)), df, int(ok.sum()), group, x, y, xs, ys, wk, ok, infl)


def weighted_corr(data: pd.DataFrame, x: str, y: str, w: np.ndarray, domain: np.ndarray,
                  how: str, method: str = "pearson", adjust: tuple[str, ...] = (),
                  label: str = "weighted", group: str = "All") -> Result:
    """Weighted correlation over a domain, with the design-based SE (rules 1-4 above).

    `data` holds every row of the pooled surveys; `domain` marks the analysed ones. Rows
    outside the domain, without the variables or with zero weight get zero influence but
    their PSUs and strata stay in the variance.
    """
    return _corr(data, x, y, w, domain, method, list(adjust), data[STRATA].to_numpy(),
                 data[PSU].to_numpy(), label, f"{how}, strata + PSU", group)


def weights_only_corr(data: pd.DataFrame, x: str, y: str, w: np.ndarray, domain: np.ndarray,
                      how: str, group: str = "All") -> Result:
    """Half-way between `weighted_corr` and `naive_corr`: the weights, but not the design.

    Same weights, so the same r as `weighted_corr`. The SE treats every row as its own PSU
    in a single stratum (rules 1-2 followed, rule 3 ignored), and df = n - 1. This is what
    a generic "weighted correlation" routine gives when it knows nothing about PSUs.
    """
    n = len(data)
    res = _corr(data, x, y, w, domain, "pearson", [], np.zeros(n), np.arange(n), "weights only",
                f"{how}, rows independent", group)
    res.df = res.n - 1
    return res


def naive_corr(data: pd.DataFrame, x: str, y: str, domain: np.ndarray, method: str = "pearson",
               adjust: tuple[str, ...] = (), group: str = "All") -> Result:
    """The naive twin: the textbook correlation of the sample, as scipy.stats gives it.

    The same estimator with every weight equal to 1 (so r is the ordinary Pearson, Spearman
    or partial correlation), then the textbook formulas for independent rows:
    SE = sqrt((1 - r^2) / df) with df = n - 2 - covariates, p from t = r / SE, and the
    Fisher-z interval tanh(atanh(r) +- 1.96 / sqrt(df - 1)). It describes the people who
    took part, not the US population, and ignores the clustering.
    """
    n = len(data)
    res = _corr(data, x, y, np.ones(n), domain, method, list(adjust), np.zeros(n),
                np.arange(n), "naive", "no weights, rows independent", group)
    res.df = res.n - 2 - len(adjust)
    res.se = float(np.sqrt((1 - res.r ** 2) / res.df))
    half = stats.norm.ppf(0.5 + CI / 2) / np.sqrt(res.df - 1)
    res.textbook_bounds = (float(np.tanh(np.arctanh(res.r) - half)), float(np.tanh(np.arctanh(res.r) + half)))
    return res


def difference(a: Result, b: Result, data: pd.DataFrame) -> Result:
    """r of group a minus r of group b, two groups of one design (men and women).

    The groups share PSUs, so their r are not independent: the influence rows are subtracted
    first and then summed per PSU, which keeps the covariance between the two estimates.
    """
    var, df = design_variance(a.infl - b.infl, data[STRATA].to_numpy(), data[PSU].to_numpy())
    return Result(a.label, a.how, a.r - b.r, float(np.sqrt(var)), df, a.n + b.n,
                  f"{a.group} minus {b.group}", is_difference=True)


def difference_independent(a: Result, b: Result) -> Result:
    """r of a minus r of b for independent samples: variances and degrees of freedom add.

    Right for two surveys (separate samples, separate strata) and what the naive analysis
    assumes for any two groups. Computed on the r scale (not Fisher's z test).
    """
    return Result(a.label, a.how, a.r - b.r, float(np.hypot(a.se, b.se)), a.df + b.df, a.n + b.n,
                  f"{a.group} minus {b.group}", is_difference=True)


# ---- output ----------------------------------------------------------------------------

RESULTS: list[dict] = []
STYLE = {"weighted": ("C0", "-"), "naive": ("C3", "-"), "mistake": ("C1", "--"),  # colour, dash
         "weights only": ("C2", ":"), "rescaled": ("C4", "--")}


def _line(res: Result) -> str:
    lo, hi = res.bounds
    name = "difference in r" if res.is_difference else "r"
    return (f"{name} = {res.r:+.3f} (95% CI {lo:+.3f} to {hi:+.3f}), SE = {res.se:.4f}, "
            f"p = {res.p:.2g}, n = {res.n:,}, df = {res.df:,}")


def show(case: int, name: str, title: str, surveys: list[str], panels: list[list[Result]],
         labels: tuple[str, str] | None = None, note: str = "") -> None:
    """Print a case, add its rows to the results table and save its figure.

    `panels` is one list of results per plot panel; the first result of a panel gives the
    points and the axis names (unless `labels` is given), every result that is not a
    difference gives a fitted line.
    """
    print(f"\nCase {case}: {title}  [{', '.join(surveys)}]")
    for res in (r for panel in panels for r in panel):
        print(f"  {res.label:<12} {res.group}: {res.how}\n      {_line(res)}")
        lo, hi = res.bounds
        RESULTS.append({"case": case, "name": name, "surveys": " ".join(surveys), "group": res.group,
                        "analysis": res.label, "how": res.how, "n": res.n, "r": res.r, "se": res.se,
                        "df": res.df, "ci_low": lo, "ci_high": hi, "p": res.p})
    if note:
        print(f"  Note: {note}")

    drawn = [p for p in panels if not p[0].is_difference]
    extra = [r for p in panels if p[0].is_difference for r in p]
    fig, axes = plt.subplots(1, len(drawn), figsize=(7.5 * len(drawn), 7.6), squeeze=False)
    rng = np.random.default_rng(0)
    for ax, panel in zip(axes[0], drawn):
        first = panel[0]
        rows = np.flatnonzero(first.ok)
        rows = rng.choice(rows, min(3000, len(rows)), replace=False)
        ax.scatter(first.xs[rows], first.ys[rows], s=5, alpha=0.25, color="0.45", linewidths=0)
        grid = np.array([first.xs[first.ok].min(), first.xs[first.ok].max()])
        for res in panel:
            colour, dash = STYLE[res.label]
            x, y, w = res.xs[res.ok], res.ys[res.ok], res.w[res.ok]
            mx, my = np.average(x, weights=w), np.average(y, weights=w)
            slope = np.average((x - mx) * (y - my), weights=w) / np.average((x - mx) ** 2, weights=w)
            ax.plot(grid, my + slope * (grid - mx), color=colour, ls=dash, lw=2, label=res.label)
        ax.set_xlabel(labels[0] if labels else first.x)
        ax.set_ylabel(labels[1] if labels else first.y)
        ax.set_title(first.group if len(drawn) > 1 else "")
        ax.legend(loc="upper left", frameon=False)
        text = "\n".join(f"{r.label}: {r.how}\n  {_line(r)}" for r in panel)
        ax.text(0, -0.13, text, transform=ax.transAxes, va="top", family="monospace", fontsize=8.5)
    foot = "\n".join(f"{r.label}, {r.group}: {_line(r)}" for r in extra)
    if note:
        foot += ("\n" if foot else "") + f"Note: {note}"
    fig.text(0.02, 0.012, foot, family="monospace", fontsize=8.5, va="bottom")
    fig.suptitle(f"Case {case}: {title}\n{', '.join(surveys)}; points: up to 3,000 analysed rows",
                 fontsize=11)
    fig.subplots_adjust(left=0.09 if len(drawn) == 1 else 0.06, right=0.98, top=0.9, bottom=0.3,
                        wspace=0.18)
    OUTPUT_DIR.mkdir(exist_ok=True)
    fig.savefig(OUTPUT_DIR / f"case_{case:02d}_{name}.png", dpi=130)
    plt.close(fig)


# ---- the cases -------------------------------------------------------------------------

def case_01_one_survey() -> None:
    """The base case: one survey, the right weight, the design-based SE.

    BMI comes from the exam, glucose from the fasting subsample, so the fasting weight is
    used: it is the only one built for the people who have both values. Every later case
    changes one thing relative to this one.
    """
    surveys = ["2017-2018"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y])
    show(1, "one_survey", "one survey", surveys,
         [[weighted_corr(data, X, Y, w, adults(data), how), naive_corr(data, X, Y, adults(data))]])


def case_02_wrong_weight() -> None:
    """Differs from case 1 by a mistake: the exam weight for a fasting-subsample variable.

    The exam weight is the natural pick for BMI and is positive for everyone examined. It
    does not correct for who was assigned to the morning session and came fasting, so it
    represents the wrong population. It also keeps people who have a glucose value but a
    zero fasting weight, so n grows.
    """
    surveys = ["2017-2018"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y])
    w_exam, how_exam = analysis_weight(data, [X, Y], force=EXAM)
    show(2, "wrong_weight", "the wrong weight (exam instead of fasting)", surveys,
         [[weighted_corr(data, X, Y, w, adults(data), how),
           weighted_corr(data, X, Y, w_exam, adults(data), how_exam, label="mistake"),
           naive_corr(data, X, Y, adults(data))]])


def case_03_weight_levels() -> None:
    """Differs from case 1 by the variables: three pairs, one for each weight level.

    Interview > exam > fasting: each level is a subsample of the previous one and has its
    own weight. The pair decides the weight (rule 1), not the researcher:
    age x sedentary time, both from the interview: interview weight;
    age x BMI, interview and exam: exam weight;
    BMI x glucose, exam and fasting subsample: fasting weight.
    Nothing is forced here: `pick_weight` reads each choice from the weights dictionary.
    """
    surveys = ["2017-2018"]
    data = load(surveys)
    panels = []
    for x, y, group in [(AGE, SEDENTARY, "interview x interview"), (AGE, X, "interview x exam"),
                        (X, Y, "exam x fasting")]:
        w, how = analysis_weight(data, [x, y])
        panels.append([weighted_corr(data, x, y, w, adults(data), how, group=group),
                       naive_corr(data, x, y, adults(data), group=group)])
    show(3, "weight_levels", "three weight levels, three pairs of variables", surveys, panels)


def case_04_weights_only() -> None:
    """Differs from case 1 by a half-way analysis: the weights are used, strata and PSUs are not.

    r equals the correct one, because r depends on the weights only. The SE is too small
    and df is n - 1 instead of 15: people of one PSU resemble each other, so 2,242 rows
    carry less information than 2,242 independent people. Weights fix the estimate;
    strata and PSUs fix its uncertainty. Both are needed.
    """
    surveys = ["2017-2018"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y])
    show(4, "weights_only", "weights without strata and PSUs", surveys,
         [[weighted_corr(data, X, Y, w, adults(data), how),
           weights_only_corr(data, X, Y, w, adults(data), how),
           naive_corr(data, X, Y, adults(data))]])


def case_05_rescaled() -> None:
    """Differs from case 1 by rescaling the weights to a mean of 1 (they sum to n, not to the population).

    A common habit ("normalised weights"). For a correlation it is harmless: a constant
    factor cancels in r and in its design-based SE, so every number equals case 1. It is
    not harmless for totals, and it does not turn weights into counts: df stays 15.
    """
    surveys = ["2017-2018"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y])
    w_mean1 = w * (w > 0).sum() / w.sum()
    correct = weighted_corr(data, X, Y, w, adults(data), how)
    rescaled = weighted_corr(data, X, Y, w_mean1, adults(data), how + " rescaled to mean 1",
                             label="rescaled")
    show(5, "rescaled", "weights rescaled to a mean of 1", surveys,
         [[correct, rescaled, naive_corr(data, X, Y, adults(data))]],
         note=f"rescaled minus weighted: r {rescaled.r - correct.r:.1e}, SE {rescaled.se - correct.se:.1e}")


def case_06_two_surveys() -> None:
    """Differs from case 1 by pooling two adjacent surveys into a 4-year estimate.

    Each survey keeps its own 2-year fasting weight, divided by 2. Strata numbers are
    unique across surveys, so the PSUs of both surveys enter the variance and df doubles.
    """
    surveys = ["2015-2016", "2017-2018"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y])
    show(6, "two_surveys", "two surveys pooled", surveys,
         [[weighted_corr(data, X, Y, w, adults(data), how), naive_corr(data, X, Y, adults(data))]])


def case_07_many_surveys() -> None:
    """Differs from case 6 by the number of surveys: seven, 2005-2018, weight / 7.

    The naive SE keeps shrinking with n; the design-based SE is limited by the number of
    PSUs, and its df is 109, not about 15,000.
    """
    surveys = [f"{y}-{y + 1}" for y in range(2005, 2018, 2)]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y])
    show(7, "many_surveys", "seven surveys pooled", surveys,
         [[weighted_corr(data, X, Y, w, adults(data), how), naive_corr(data, X, Y, adults(data))]])


def case_08_pre_post_covid() -> None:
    """Differs from case 1 by comparing two surveys instead of pooling them.

    2017-2018 (before the pandemic) and 2021-2023 (after) are analysed separately, each
    with its own undivided fasting weight and its own design (df = 15 each). The surveys
    are independent samples, so the variance of the difference in r is the sum of the two
    variances and its df is 30.
    """
    surveys = ["2017-2018", "2021-2023"]
    wtd, nai = [], []
    for s in surveys:
        data = load([s])
        w, how = analysis_weight(data, [X, Y])
        wtd.append(weighted_corr(data, X, Y, w, adults(data), how, group=s))
        nai.append(naive_corr(data, X, Y, adults(data), group=s))
    show(8, "pre_post_covid", "before and after the pandemic, compared", surveys,
         [[wtd[0], nai[0]], [wtd[1], nai[1]],
          [difference_independent(*wtd), difference_independent(*nai)]])


def case_09_across_covid() -> None:
    """Differs from case 6 by pooling two surveys that are not adjacent and differ in design.

    The arithmetic is the same (own weight / 2), but 2021-2023 was sampled and weighted
    differently after the pandemic and the years 2019-2020 are missing, so the result is an
    average of two periods, not a 2017-2023 estimate. The 2017-March 2020 pre-pandemic file
    is not in this snapshot; it would replace 2017-2018, never be added to it.
    """
    surveys = ["2017-2018", "2021-2023"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y])
    show(9, "across_covid", "pooled across the pandemic", surveys,
         [[weighted_corr(data, X, Y, w, adults(data), how), naive_corr(data, X, Y, adults(data))]],
         note="two periods and two sample designs are averaged; compare with case 8")


def case_10_weight_changes() -> None:
    """Differs from case 9 by the variable: BMI x HbA1c, whose weight changes between surveys.

    HbA1c uses the exam weight in 2017-2018 and the new blood-draw weight in 2021-2023
    (introduced because many examined people gave no blood). The weight must be picked per
    survey; one column for both surveys is a mistake (here: the exam weight everywhere).
    """
    surveys = ["2017-2018", "2021-2023"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, HBA1C])
    w_exam, how_exam = analysis_weight(data, [X, HBA1C], force=EXAM)
    show(10, "weight_changes", "the weight changes between surveys (BMI x HbA1c)", surveys,
         [[weighted_corr(data, X, HBA1C, w, adults(data), how),
           weighted_corr(data, X, HBA1C, w_exam, adults(data), how_exam, label="mistake"),
           naive_corr(data, X, HBA1C, adults(data))]])


def _subgroup_case(case: int, name: str, title: str, group: str, member) -> None:
    """Cases 11 and 12: one subgroup as a domain (correct) and with the other rows deleted first."""
    surveys = ["2017-2018"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y])
    subgroup = adults(data) & member(data).to_numpy()
    kept = data[subgroup].reset_index(drop=True)
    deleted_first = weighted_corr(kept, X, Y, w[subgroup], np.ones(len(kept), bool),
                                  how + ", other rows deleted", label="mistake", group=group)
    show(case, name, title, surveys,
         [[weighted_corr(data, X, Y, w, subgroup, how, group=group), deleted_first,
           naive_corr(data, X, Y, subgroup, group=group)]])


def case_11_subgroup() -> None:
    """Differs from case 1 by a filter: Mexican American adults, analysed as a domain.

    Correct: all rows of the survey stay and the subgroup is a mask (rule 4). Mistake:
    delete the other rows first (`data[data.Ethnicity == 1]`). r is the same, but PSUs with
    no Mexican American participant vanish and strata left with one PSU add no variance,
    so the SE shrinks and df falls from 15 to 11.
    """
    _subgroup_case(11, "subgroup", "a subgroup (Mexican American adults)", "Mexican American",
                   lambda data: data[ETHNICITY] == 1)


def case_12_age_band() -> None:
    """Differs from case 11 by the filter: age 40-59, a group found in every PSU.

    The same two analyses. Here deleting the rows first loses no PSU and no stratum, so
    the "mistake" gives the same SE and df as the domain analysis. That is luck, not
    correctness: whether deletion is harmless depends on how the group is spread over the
    PSUs, and the domain analysis is right in both cases.
    """
    _subgroup_case(12, "age_band", "a subgroup (age 40-59)", "Age 40-59",
                   lambda data: data[AGE].between(40, 59))


def case_13_by_gender() -> None:
    """Differs from case 11 by having two domains of one design and comparing them.

    Men and women live in the same PSUs, so their two r are not independent: `difference`
    takes the SE from the joint design (df = 15). The naive twin treats them as two
    independent samples.
    """
    surveys = ["2017-2018"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y])
    groups = {"Men": adults(data) & (data[GENDER] == 1).to_numpy(),
              "Women": adults(data) & (data[GENDER] == 0).to_numpy()}
    wtd = [weighted_corr(data, X, Y, w, d, how, group=g) for g, d in groups.items()]
    nai = [naive_corr(data, X, Y, d, group=g) for g, d in groups.items()]
    show(13, "by_gender", "by gender, and the difference", surveys,
         [[wtd[0], nai[0]], [wtd[1], nai[1]],
          [difference(*wtd, data), difference_independent(*nai)]])


def case_14_adjusted() -> None:
    """Differs from case 1 by adjusting for age and gender (partial correlation).

    BMI and glucose are each regressed on age and gender by weighted least squares, and the
    residuals are correlated with the same weight and design. The naive twin uses ordinary
    least squares. The SE treats the residuals as fixed.
    """
    surveys = ["2017-2018"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y, AGE, GENDER])
    adjust = (AGE, GENDER)
    show(14, "adjusted", "adjusted for age and gender", surveys,
         [[weighted_corr(data, X, Y, w, adults(data), how, adjust=adjust),
           naive_corr(data, X, Y, adults(data), adjust=adjust)]],
         labels=(f"{X}, residual", f"{Y}, residual"))


def case_15_spearman() -> None:
    """Differs from case 1 by the method: Spearman, because glucose is strongly skewed.

    Ranks are weighted too: each value is replaced by the share of the population below it
    (plus half of its ties), then correlated with the same weight and design. The naive
    twin is the ordinary Spearman correlation of sample ranks.
    """
    surveys = ["2017-2018"]
    data = load(surveys)
    w, how = analysis_weight(data, [X, Y])
    show(15, "spearman", "Spearman (rank) correlation", surveys,
         [[weighted_corr(data, X, Y, w, adults(data), how, method="spearman"),
           naive_corr(data, X, Y, adults(data), method="spearman")]],
         labels=(f"{X}, rank (0-1)", f"{Y}, rank (0-1)"))


def case_16_refused() -> None:
    """Differs from all others: two analyses that must be refused, not computed.

    1999-2000 + 2001-2002 need the 4-year weights NCHS issued for that pair; dividing the
    2-year weights by 2 is wrong there. 2017-2020 (pre-pandemic) is not in the snapshot.
    A refusal is the correct result; an r would be a defect.
    """
    print("\nCase 16: analyses that must be refused")
    attempts = {"1999-2000 2001-2002": lambda: analysis_weight(load(["1999-2000", "2001-2002"]), [X, Y]),
                "2017-2020": lambda: load(["2017-2020"])}
    for surveys, attempt in attempts.items():
        try:
            attempt()
            outcome = "NOT REFUSED (defect)"
        except Refused as why:
            outcome = f"refused: {why}"
        print(f"  {surveys}: {outcome}")
        RESULTS.append({"case": 16, "name": "refused", "surveys": surveys, "analysis": outcome})


CASES = [case_01_one_survey, case_02_wrong_weight, case_03_weight_levels, case_04_weights_only,
         case_05_rescaled, case_06_two_surveys, case_07_many_surveys, case_08_pre_post_covid,
         case_09_across_covid, case_10_weight_changes, case_11_subgroup, case_12_age_band,
         case_13_by_gender, case_14_adjusted, case_15_spearman, case_16_refused]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="*", type=int, help="case numbers (default: all)")
    ap.add_argument("--list", action="store_true", help="print the case list and stop")
    args = ap.parse_args()
    if args.list:
        for i, case in enumerate(CASES, 1):
            print(f"{i:>2}  {case.__name__:<26} {case.__doc__.splitlines()[0]}")
        return
    bad = fetch_data.problems()
    if bad:
        sys.exit("data/ is not ready (" + "; ".join(bad) + "). Run: uv run python fetch_data.py")
    for i in args.cases or range(1, len(CASES) + 1):
        CASES[i - 1]()
    OUTPUT_DIR.mkdir(exist_ok=True)
    fields = ["case", "name", "surveys", "group", "analysis", "how", "n", "r", "se", "df",
              "ci_low", "ci_high", "p"]
    with open(OUTPUT_DIR / "results.csv", "w", newline="") as f:
        out = csv.DictWriter(f, fieldnames=fields)
        out.writeheader()
        out.writerows(RESULTS)
    print(f"\nFigures and results.csv are in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
