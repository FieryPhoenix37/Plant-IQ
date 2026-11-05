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
    pairwise,
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


def jitter(ax, lf: pl.LazyFrame, by: str, order: list, labels: dict | None = None,
           outcome: str = "school", seed: int = 0,
           ylabel: str = "School transmission cases\n(% of school population)") -> None:
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

    ax.axhline(SIZABLE_OUTBREAK * 100, color="gray", ls="--", lw=1)
    # % of runs > 5% goes in the tick label (as printed under each column in the paper)
    ticks = [f"{(labels or {}).get(c, c)}\n{stats[c]['pct_gt5']:.1f}%" for c in order]
    ax.set_xticks(range(len(order)), ticks)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylabel(ylabel)


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
        ax.set_xlabel("Introduction rate")
    # Only after every panel is drawn: setting limits earlier freezes the shared y-axis.
    np.atleast_1d(axes)[0].set_ylim(bottom=-1)
    return _save(fig, f"fig3_jitter_{school.name}.png")


def fig5_jitter(school: School) -> Path:
    lf = figure_frame(school, "fig5_jitter")
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    jitter(ax, lf, "cohorting", COHORTING_ORDER, COHORTING_LABELS)
    ax.set_title(f"{school.name.title()} school, R0 = 1.5, weekly introductions, no testing")
    ax.set_xlabel("Cohorting")
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


def fig6_pairwise(school: School) -> Path:
    """Pairwise ln(mean outbreak size) ratios, column vs row, with Mann-Whitney U significance."""
    keys = [(c, t) for c in COHORTING_ORDER for t in TESTING_ORDER]
    pw = pairwise(figure_frame(school, "fig6_pairwise"), ["cohorting", "testing"], "school")
    cell = {((r["row_cohorting"], r["row_testing"]), (r["col_cohorting"], r["col_testing"])): r
            for r in pw.iter_rows(named=True)}
    grid = np.array([[cell[(rk, ck)]["log_ratio"] for ck in keys] for rk in keys])
    lim = np.nanmax(np.abs(grid))

    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(grid, cmap="RdBu_r", vmin=-lim, vmax=lim)
    for i, rk in enumerate(keys):
        for j, ck in enumerate(keys):
            p = cell[(rk, ck)]["p_value"]
            mark = "**" if p < 0.01 else "*" if p < 0.05 else ""
            if mark and i != j:
                ax.text(j, i, mark, ha="center", va="center", fontsize=7)
    label = [f"{COHORTING_LABELS[c]} / {TESTING_LABELS[t]}" for c, t in keys]
    ax.set_xticks(range(len(keys)), label, rotation=90, fontsize=7)
    ax.set_yticks(range(len(keys)), label, fontsize=7)
    ax.set_title(f"{school.name.title()}: ln(mean outbreak size, column / row)\n"
                 "R0 = 1.5, weekly introductions; ** p < 0.01, * p < 0.05 (Mann-Whitney U)", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.8, label="blue: column better, red: column worse")
    return _save(fig, f"fig6_pairwise_{school.name}.png")


def fig7_jitter(school: School) -> Path:
    """No vaccination vs all teachers/staff vaccinated, by cohorting, students and teachers+staff."""
    lf = figure_frame(school, "fig7_jitter").with_columns(
        group=pl.format("{}|{}", pl.col("cohorting"), pl.col("teacher_vaccine").cast(pl.Int8))
    )
    order = [f"{c}|{v}" for c in COHORTING_ORDER for v in (0, 1)]
    labels = {f"{c}|{v}": f"{COHORTING_LABELS[c]}\n{'vaccinated' if v else 'unvaccinated'}"
              for c in COHORTING_ORDER for v in (0, 1)}
    fig, axes = plt.subplots(2, 1, figsize=(10, 9.5), constrained_layout=True)
    for ax, outcome, who in zip(axes, ["students", "teachers_staff"], ["students", "teachers + staff"]):
        jitter(ax, lf, "group", order, labels, outcome=outcome,
               ylabel=f"School transmission cases\n(% of {who})")
        ax.tick_params(axis="x", labelsize=8)
        ax.set_ylim(bottom=-1)
        ax.set_title(f"{school.name.title()}: {who} (R0 = 1.5, weekly introductions, no testing)", fontsize=10)
    return _save(fig, f"fig7_jitter_{school.name}.png")


def fig2_network(school: School, seed: int = 1) -> Path:
    """Example contact network on a day when only student cohort 1 is on campus (Fig. 2C/D).

    Illustrative: uses the upstream generators with their default settings, which are not yet fitted
    to the Zenodo network statistics (Part 1B).
    """
    import networkx as nx

    from .seirsplus_k12.networks import generate_highschool_contact_network, generate_K5_school_contact_network

    np.random.seed(seed)  # upstream generators use the global numpy RNG
    if school.name == "primary":
        nets, info = generate_K5_school_contact_network(
            num_grades=6, num_classrooms_per_grade=4, class_sizes=20, num_student_blocks=2, num_staff=24)
    else:
        nets, info = generate_highschool_contact_network(num_student_blocks=2)
    G = nets["onsite-block1"].copy()
    # The K-5 generator leaves weight-0 self-loops (stored zeros from scipy block_diag); not contacts.
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    pos = nx.spring_layout(G, weight="layout_weight", seed=seed, iterations=60)

    role = np.array(info["node_labels"])
    grade = {s: g for g, ids in info["grades_studentIDs"].items() for s in ids}
    offsite = {s for s, b in info["studentIDs_studentBlocks"].items() if b != 1}
    blues = plt.get_cmap("Blues")
    colors = []
    for n in G.nodes():
        if role[n] == "student":
            colors.append(blues(0.45 + 0.5 * grade[n] / max(grade.values())))
        else:
            colors.append("tab:green" if role[n] == "teacher" else "gold")

    fig, ax = plt.subplots(figsize=(7, 7))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="lightgray", alpha=0.25, width=0.3)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=[6 if n in offsite else 10 for n in G.nodes()],
                           node_color=colors, linewidths=0)
    ax.set_title(f"{school.name.title()} school, cohort 1 on campus\n"
                 "blue: students (shade = grade), green: teachers, yellow: staff;\n"
                 "peripheral nodes: off-campus cohort (household edges only)", fontsize=9)
    ax.set_axis_off()
    return _save(fig, f"fig2_network_{school.name}.png")


def make_all(school: School) -> list[Path]:
    figs = (fig2_network, fig3_jitter, fig3_heatmap, fig4_heatmap, fig5_jitter, fig6_pairwise, fig7_jitter)
    return [f(school) for f in figs]


def _save(fig, name: str) -> Path:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path
