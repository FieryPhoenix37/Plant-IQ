"""Part 1A: evidence behind the data-handling decisions in zenodo.py, written to results/.

- zenodo_anomalies.csv: per school, counts of each known data anomaly (see zenodo.anomalies).
- zenodo_block_comparison.csv: for every secondary scenario with two 1,000-run blocks, the
  canonical (high-clustering network) vs alt_network block side by side: network statistics,
  outbreak sizes, share of runs > 5% with a binomial z-score for the difference, and counts above
  the population (all of which lie in alt_network blocks).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .zenodo import RESULTS_DIR, SCENARIO_COLS, SIZABLE_OUTBREAK, School, anomalies, load


def anomaly_table(schools: list[School]) -> pl.DataFrame:
    return pl.DataFrame([{"school": s.name, **anomalies(s)} for s in schools])


def block_comparison(school: School) -> pl.DataFrame:
    """One row per two-block scenario: canonical vs alt_network statistics and differences."""
    keys = list(SCENARIO_COLS.values())
    over = (pl.col("n_infected_students") > school.n_students) | (
        pl.col("n_infected_teachers_staff") > school.n_teachers_staff
    )
    per_block = (
        load(school)
        .filter(pl.col("block").n_unique().over(keys) == 2)
        .group_by(*keys, "version")
        .agg(
            clustering=pl.col("clustering").mean(),
            deg_mean=pl.col("deg_mean").mean(),
            mean_infected=pl.col("n_infected").mean(),
            mean_peak=pl.col("n_peak_infected").mean(),
            median_frac=pl.col("frac_endo").median(),
            pct_gt5=(pl.col("frac_endo") > SIZABLE_OUTBREAK).mean() * 100,
            n_above_population=over.sum(),
        )
        .collect()
    )
    if per_block.is_empty():
        return per_block
    stats = ["clustering", "deg_mean", "mean_infected", "mean_peak", "median_frac", "pct_gt5", "n_above_population"]
    # pivot orders new columns by first appearance, which parallel group_by does not fix: pin it.
    wide = per_block.pivot(on="version", index=keys, values=stats).select(
        *keys, *(f"{s}_{v}" for s in stats for v in ("canonical", "alt_network"))
    )
    p = (pl.col("pct_gt5_alt_network") + pl.col("pct_gt5_canonical")) / 200
    se = (p * (1 - p) * 2 / 1000).sqrt() * 100  # SE of a difference of two 1,000-run proportions
    return wide.with_columns(
        infected_ratio=pl.col("mean_infected_alt_network") / pl.col("mean_infected_canonical"),
        pct_gt5_diff=pl.col("pct_gt5_alt_network") - pl.col("pct_gt5_canonical"),
    ).with_columns(
        pct_gt5_z=pl.when(se > 0).then(pl.col("pct_gt5_diff") / se).otherwise(0.0),
    ).sort(keys)


def write_all(schools: list[School]) -> list[Path]:
    RESULTS_DIR.mkdir(exist_ok=True)
    out = [RESULTS_DIR / "zenodo_anomalies.csv"]
    anomaly_table(schools).write_csv(out[0])
    blocks = [block_comparison(s).with_columns(school=pl.lit(s.name)) for s in schools]
    blocks = [b for b in blocks if b.height]
    if blocks:
        out.append(RESULTS_DIR / "zenodo_block_comparison.csv")
        pl.concat(blocks, how="diagonal_relaxed").write_csv(out[-1], float_precision=4)
    return out
