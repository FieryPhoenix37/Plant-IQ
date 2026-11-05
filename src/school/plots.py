"""Part 1A figures in the style of McGee et al. 2021: jitter plots (Figs. 3A-D, 5) and heatmaps (Figs. 3E-F, 4)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from .zenodo import (
    COHORTING_ORDER,
    INTRO_ORDER,
    OUTCOMES,
    ROOT,
    SIZABLE_OUTBREAK,
    TESTING_ORDER,
    School,
    figure_frame,
    summarize,
)

FIG_DIR = ROOT / "figures" / "part1a"
MEDIAN_COLOR, P95_COLOR = "black", "tab:orange"

COHORTING_LABELS = {
    "all-onsite": "No cohorting",
    "alternate-blocks-daily": "Alternating days",
    "alternate-blocks-weekly": "Alternating weeks",
}
TESTING_LABELS = {
    "none": "No testing",
    "weekly_teachers-only": "Weekly T",
    "semiweekly_teachers-only": "2x/wk T",
    "weekly_teachers-and-students": "Weekly T+S",
    "semiweekly_teachers-and-students": "2x/wk T+S",
}


def jitter(ax, lf: pl.LazyFrame, by: str, order: list[str], labels: dict | None = None,
           outcome: str = "school", seed: int = 0) -> None:
    """One jitter column per category: points, median (black), 95th pct (orange), % runs > 5% below."""
    rng = np.random.default_rng(seed)
    col = OUTCOMES[outcome]
    stats = {r[by]: r for r in summarize(lf, [by], outcome).iter_rows(named=True)}
    data = lf.select(by, col).collect()

    for i, cat in enumerate(order):
        y = data.filter(pl.col(by) == cat)[col].to_numpy() * 100
        x = i + rng.uniform(-0.3, 0.3, y.size)
        ax.scatter(x, y, s=4, alpha=0.35, color="tab:blue", linewidths=0)
        s = stats[cat]
        ax.hlines(s["median"] * 100, i - 0.35, i + 0.35, color=MEDIAN_COLOR, lw=2)
        ax.hlines(s["p95"] * 100, i - 0.35, i + 0.35, color=P95_COLOR, lw=2)
        ax.annotate(f"{s['pct_gt5']:.1f}%", (i, 0), xytext=(0, -28), textcoords="offset points",
                    xycoords=("data", "axes fraction"), ha="center", fontsize=9)

    ax.axhline(SIZABLE_OUTBREAK * 100, color="gray", ls="--", lw=1)
    ax.set_xticks(range(len(order)), [(labels or {}).get(c, c) for c in order])
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylabel("School transmission cases\n(% of school population)")


def heatmap(ax, table: pl.DataFrame, rows: list[str], cols: list[str], row_order, col_order,
            title: str, row_fmt=str, col_fmt=str):
    """Fraction of runs > 5% as a heatmap. row_order / col_order are lists of key tuples."""
    lookup = {(tuple(r[k] for k in rows), tuple(r[k] for k in cols)): r["pct_gt5"] / 100
              for r in table.iter_rows(named=True)}
    grid = np.array([[lookup.get((ro, co), np.nan) for co in col_order] for ro in row_order])
    im = ax.imshow(grid, vmin=0, vmax=1, cmap="magma_r", aspect="auto")
    ax.set_yticks(range(len(row_order)), [row_fmt(r) for r in row_order], fontsize=7)
    ax.set_xticks(range(len(col_order)), [col_fmt(c) for c in col_order], rotation=90, fontsize=7)
    ax.set_title(title, fontsize=10)
    return im


def _r0_values(lf: pl.LazyFrame) -> list[float]:
    return sorted(lf.select(pl.col("R0").unique()).collect()["R0"].to_list())


def fig3_jitter(school: School) -> Path:
    lf = figure_frame(school, "fig3_jitter")
    r0s = _r0_values(lf)
    fig, axes = plt.subplots(1, len(r0s), figsize=(5 * len(r0s), 4.5), sharey=True)
    for ax, r0 in zip(np.atleast_1d(axes), r0s):
        jitter(ax, lf.filter(pl.col("R0") == r0), "intro", INTRO_ORDER)
        ax.set_title(f"{school.name.title()} school, R0 = {r0}")
        ax.set_xlabel("Introduction rate", labelpad=22)
    # Only after every panel is drawn: setting limits earlier freezes the shared y-axis.
    np.atleast_1d(axes)[0].set_ylim(bottom=-1)
    return _save(fig, f"fig3_jitter_{school.name}.png")


def fig5_jitter(school: School) -> Path:
    lf = figure_frame(school, "fig5_jitter")
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    jitter(ax, lf, "cohorting", COHORTING_ORDER, COHORTING_LABELS)
    ax.set_title(f"{school.name.title()} school, R0 = 1.5, weekly introductions, no testing")
    ax.set_xlabel("Cohorting", labelpad=22)
    ax.set_ylim(bottom=-1)
    return _save(fig, f"fig5_jitter_{school.name}.png")


def fig3_heatmap(school: School) -> Path:
    lf = figure_frame(school, "fig3_heatmap")
    r0s = _r0_values(lf)
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True)
    for ax, outcome in zip(axes, ["students", "teachers_staff"]):
        im = heatmap(ax, summarize(lf, ["R0", "intro"], outcome), ["R0"], ["intro"],
                     [(r,) for r in r0s], [(i,) for i in INTRO_ORDER],
                     f"{school.name.title()}: {outcome.replace('_', ' + ')}",
                     row_fmt=lambda r: f"R0={r[0]}", col_fmt=lambda c: c[0])
    fig.colorbar(im, ax=axes, label="Fraction of runs > 5% infected")
    return _save(fig, f"fig3_heatmap_{school.name}.png")


def fig4_heatmap(school: School) -> Path:
    lf = figure_frame(school, "fig4_heatmap")
    keys_r, keys_c = ["R0", "intro"], ["cohorting", "testing"]
    rows = [(r, i) for r in _r0_values(lf) for i in INTRO_ORDER]
    cols = [(c, t) for c in COHORTING_ORDER for t in TESTING_ORDER]
    fig, axes = plt.subplots(1, 2, figsize=(12, 0.22 * len(rows) + 2.5), constrained_layout=True)
    for ax, outcome in zip(axes, ["students", "teachers_staff"]):
        im = heatmap(ax, summarize(lf, keys_r + keys_c, outcome), keys_r, keys_c, rows, cols,
                     f"{school.name.title()}: {outcome.replace('_', ' + ')}",
                     row_fmt=lambda r: f"R0={r[0]}, {r[1]}",
                     col_fmt=lambda c: f"{COHORTING_LABELS[c[0]]} / {TESTING_LABELS[c[1]]}")
    fig.colorbar(im, ax=axes, label="Fraction of runs > 5% infected")
    return _save(fig, f"fig4_heatmap_{school.name}.png")


def make_all(school: School) -> list[Path]:
    return [f(school) for f in (fig3_jitter, fig3_heatmap, fig4_heatmap, fig5_jitter)]


def _save(fig, name: str) -> Path:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path
