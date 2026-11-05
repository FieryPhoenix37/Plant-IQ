"""Part 1A: rebuild McGee et al. 2021 (PNAS) Fig. 3-5 targets from the Zenodo output (10.5281/zenodo.5237328).

Pipeline:
    1. extract():   one streaming DuckDB pass over each multi-GB .csv.gz -> slim parquet in data/
    2. targets():   filter to the paper's baseline, summarize per figure -> results/zenodo_targets.csv
    3. school.plots draws the figures from the same parquet.

Facts about the raw files, verified against the data (not the Zenodo description):
    - Comma-separated (not tab), with leftover pandas index columns and quoted tuple fields.
    - ``*_pct_*`` outcome columns are fractions in [0, 1], so ">5%" means ``> 0.05``.
    - Per group, in every row: endogenous = infected - introductions - 2 for students and - 0 for
      teachers/staff. The initial seed case is not in ``total_num_introductions``, so one case too
      many is subtracted (single-seed runs show -1). We keep the published columns as the primary
      outcome and expose a seed-corrected variant (+1 case) for sensitivity.
    - ``total_num_infected`` can exceed students + teachers/staff by 0-5 infections, more often with
      more introductions; the school-wide endogenous count inherits that excess.
    - The secondary school-wide percentage divides by 925 (= 800 students + 125 teachers), although
      the teachers/staff percentage divides by 200 (125 + 75). We use the denominators as published.
    - Rows come in blocks of 1,000 runs (``row_idx // 1000``). 2,546 secondary scenarios have two
      blocks: same inputs and network statistics, same share of runs > 5% (within noise), but
      outbreak sizes up to ~2.3x apart, growing with outbreak size. The larger block holds every
      count above the population (up to 1,880 of 1,000). The smaller block is the complete version
      (it also covers the 334 single-block scenarios), so it is labelled canonical and used for the
      targets; see label_versions().

Baseline (paper main text): student susceptibility 0.6 (primary) / 1.0 (secondary), only positive
individuals isolated (no classroom quarantine), no vaccination.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
ZENODO_FILE_URL = "https://zenodo.org/api/records/5237328/files/{name}/content"

SIZABLE_OUTBREAK = 0.05  # strictly greater than 5% of the (sub)population
RUNS_PER_BLOCK = 1000  # the authors wrote results in blocks of 1,000 replicates of one scenario


@dataclass(frozen=True)
class School:
    name: str
    raw_file: str
    md5: str
    # Denominators of the published *_pct_* columns (verified against the counts in validate()).
    total_denominator: int
    n_students: int
    n_teachers_staff: int
    susceptibility: float  # baseline relative susceptibility of students

    @property
    def raw_path(self) -> Path:
        return DATA_DIR / self.raw_file

    @property
    def parquet_path(self) -> Path:
        return DATA_DIR / f"zenodo_{self.name}.parquet"


SCHOOLS = {
    "primary": School(
        "primary", "simData_primarySchools.csv.gz", "422e4f9d1d8e65e7f206ddf64cec2439",
        total_denominator=528, n_students=480, n_teachers_staff=48, susceptibility=0.6,
    ),
    "secondary": School(
        "secondary", "simData_secondarySchools.csv.gz", "d5612dc07e98110817bde1a3a1ad6f00",
        total_denominator=925, n_students=800, n_teachers_staff=200, susceptibility=1.0,
    ),
}

# Raw column -> short name. Only these are kept in the parquet. SCENARIO_COLS is every input
# parameter in column_descriptions.csv (all columns before the network statistics), so the scenario
# key is complete: a sample-based subset missed a vaccine parameter that varies only when
# teachers are vaccinated.
SCENARIO_COLS = {
    "R0_MEAN": "R0",
    "R0_CV": "R0_cv",
    "SUSCEPTIBILITY_students": "susceptibility",
    "AVG_INTRODUCTIONS_PER_DAY": "intro_rate",
    "ONSITE_CADENCE": "cohorting",
    "TESTING_PROTOCOL": "testing",
    "TESTING_CADENCE": "testing_cadence",
    "TESTING_COMPLIANCE_RATE_RANDOM_students": "test_compliance_students",
    "TESTING_COMPLIANCE_RATE_RANDOM_teachers": "test_compliance_teachers",
    "ISOLATION_COMPLIANCE_RATE_POSITIVE_GROUPMATE": "classroom_isolation",
    "ISOLATION_COMPLIANCE_RATE_POSITIVE_CONTACT": "contact_isolation",
    "INIT_VACCINE_UPTAKE_students": "student_vaccine",
    "INIT_VACCINE_UPTAKE_teachers": "teacher_vaccine",
    "VACCINE_EFFICACY": "vaccine_efficacy",
    "VACCINE_SUSCEPTIBILITY_FACTOR": "vaccine_susceptibility_factor",
    "VACCINE_TRANSMISSIBILITY_FACTOR": "vaccine_transmissibility_factor",
    "STOP_INTERVENTION_VACCINE_DOSE": "stop_intervention_vaccine_dose",
}
# Leftover pandas indices from the authors' concatenation of result files. "Unnamed: 0" is the global
# row number (identical to the blank-named first column, whose DuckDB name depends on the column
# count, so it is not used); "Unnamed: 0.1.1" is the run's position within its 1,000-run block.
PROVENANCE_COLS = {
    "Unnamed: 0": "row_idx",
    "Unnamed: 0.1.1": "idx_in_block",
}
# Per-run statistics of the generated contact network (targets for fitting generators in Part 1B).
NETWORK_COLS = {
    "degree_mean": "deg_mean",
    "degree_CV": "deg_cv",
    "cluster_coeff": "clustering",
    "assortativity": "assortativity",
    "degree_mean_students": "deg_mean_students",
    "cluster_coeff_students": "clustering_students",
    "degree_mean_teachers-staff": "deg_mean_teachers_staff",
}
COUNT_COLS = {
    "runID": "run_id",
    "total_num_infected": "n_infected",
    "total_num_infected_students": "n_infected_students",
    "total_num_infected_teachers-staff": "n_infected_teachers_staff",
    "total_num_introductions": "n_introductions",
    "total_num_introductions_students": "n_introductions_students",
    "total_num_introductions_teachers": "n_introductions_teachers",
    "peak_num_infected": "n_peak_infected",
    "total_num_endogenous_cases": "n_endo",
    "total_num_endogenous_cases_students": "n_endo_students",
    "total_num_endogenous_cases_teachers-staff": "n_endo_teachers_staff",
    "total_pct_endogenous_cases": "frac_endo",
    "total_pct_endogenous_cases_students": "frac_endo_students",
    "total_pct_endogenous_cases_teachers-staff": "frac_endo_teachers_staff",
}

# Display order and labels for the categorical axes used in the paper.
INTRO_LABELS = {0: "single", 30: "monthly", 7: "weekly", 1: "daily"}  # keyed by round(1/rate)
COHORTING_ORDER = ["all-onsite", "alternate-blocks-daily", "alternate-blocks-weekly"]
TESTING_ORDER = [
    "none",
    "weekly_teachers-only",
    "semiweekly_teachers-only",
    "weekly_teachers-and-students",
    "semiweekly_teachers-and-students",
]
INTRO_ORDER = list(INTRO_LABELS.values())

# Outcome name -> (column, population label). "_seedfix" variants add back the over-subtracted case.
OUTCOMES = {
    "school": "frac_endo",
    "students": "frac_endo_students",
    "teachers_staff": "frac_endo_teachers_staff",
    "school_seedfix": "frac_endo_seedfix",
    "students_seedfix": "frac_endo_students_seedfix",
}


# --------------------------------------------------------------------------------------------
# 1. Extraction
# --------------------------------------------------------------------------------------------

def extract(school: School, overwrite: bool = False) -> Path:
    """Stream the raw .csv.gz once with DuckDB and write only the needed columns to parquet."""
    import duckdb

    out = school.parquet_path
    if out.exists() and not overwrite:
        return out
    if not school.raw_path.exists():
        raise FileNotFoundError(
            f"{school.raw_path} missing; download from {ZENODO_FILE_URL.format(name=school.raw_file)}"
        )

    select = ", ".join(f'"{raw}" AS {short}' for raw, short in {**PROVENANCE_COLS, **SCENARIO_COLS, **NETWORK_COLS, **COUNT_COLS}.items())
    tmp = out.with_suffix(".parquet.tmp")
    con = duckdb.connect()
    con.execute(
        f"""
        COPY (
            SELECT {select}
            FROM read_csv('{school.raw_path}', header = true, delim = ',', quote = '"',
                          compression = 'gzip', sample_size = 100000)
        ) TO '{tmp}' (FORMAT parquet, COMPRESSION zstd)
        """
    )
    con.close()
    tmp.rename(out)
    return out


# --------------------------------------------------------------------------------------------
# 2. Loading, validation, baseline filtering
# --------------------------------------------------------------------------------------------

def load(school: School) -> pl.LazyFrame:
    """Lazy frame over the slim parquet with derived columns (block/version, intro label, seed fix)."""
    lf = pl.scan_parquet(school.parquet_path).with_columns(
        intro=_intro_label(pl.col("intro_rate")),
        frac_endo_seedfix=(pl.col("n_endo") + 1) / school.total_denominator,
        frac_endo_students_seedfix=(pl.col("n_endo_students") + 1) / school.n_students,
    )
    return label_versions(lf)


def label_versions(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Add ``block`` and ``version`` ("canonical" / "inflated").

    Every primary scenario, and 334 secondary ones, has one block of 1,000 runs. The other 2,546
    secondary scenarios have two blocks with identical inputs and networks but outbreak sizes up to
    ~2.3x apart (the larger block holds every count above the population; the share of runs > 5% is
    the same in both). The lower-mean block is labelled canonical (ties: lower block number). At low
    R0 the two blocks are statistically indistinguishable, so the choice there is immaterial.
    """
    keys = list(SCENARIO_COLS.values())
    block_mean = pl.col("n_infected").mean().over("block")
    return (
        lf.with_columns(block=pl.col("row_idx") // RUNS_PER_BLOCK)
        # Ties (2 secondary scenarios with ~1-8 mean infections) go to the lower block number.
        .with_columns(rank=pl.struct(block_mean.alias("m"), pl.col("block")).rank("dense").over(keys))
        .with_columns(version=pl.when(pl.col("rank") == 1).then(pl.lit("canonical")).otherwise(pl.lit("inflated")))
        .drop("rank")
    )


def _intro_label(rate: pl.Expr) -> pl.Expr:
    """Map stored float32 rates (0, 0.0333.., 0.1428.., 1.0) to labels via the period in days."""
    period = pl.when(rate > 0).then((1 / rate).round(0).cast(pl.Int64)).otherwise(0)
    return period.replace_strict(INTRO_LABELS, return_dtype=pl.Utf8)


def validate(school: School) -> pl.DataFrame:
    """Check the invariants the analysis relies on. Raises on violation; returns scenario sizes."""
    lf = load(school)
    # Each expression must take exactly one value across all rows.
    rules = {
        "total: infected - introductions - endogenous": (
            pl.col("n_infected") - pl.col("n_introductions") - pl.col("n_endo"), [2]),
        "students: infected - introductions - endogenous": (
            pl.col("n_infected_students") - pl.col("n_introductions_students") - pl.col("n_endo_students"), [2]),
        "teachers/staff: infected - introductions - endogenous": (
            pl.col("n_infected_teachers_staff") - pl.col("n_introductions_teachers")
            - pl.col("n_endo_teachers_staff"), [0]),
        "introductions: total - students - teachers/staff": (
            pl.col("n_introductions") - pl.col("n_introductions_students") - pl.col("n_introductions_teachers"), [0]),
    }
    for label, (expr, expected) in rules.items():
        got = lf.select(expr.unique().alias("v")).collect()["v"].to_list()
        if got != expected:
            raise AssertionError(f"{school.name}: {label} = {got}, expected {expected}")
    excess = pl.col("n_infected") - pl.col("n_infected_students") - pl.col("n_infected_teachers_staff")
    if lf.select(excess.min()).collect().item() < 0:
        raise AssertionError(f"{school.name}: total infected is below the sum of its groups")

    frac_ok = lf.select(
        (pl.col("n_endo") / school.total_denominator - pl.col("frac_endo")).abs().max().alias("total"),
        (pl.col("n_endo_students") / school.n_students - pl.col("frac_endo_students")).abs().max().alias("stud"),
        (pl.col("n_endo_teachers_staff") / school.n_teachers_staff
         - pl.col("frac_endo_teachers_staff")).abs().max().alias("ts"),
    ).collect()
    if frac_ok.max_horizontal().item() > 1e-9:
        raise AssertionError(f"{school.name}: population denominators don't match: {frac_ok}")

    keys = list(SCENARIO_COLS.values())
    blocks = lf.group_by("block").agg(
        n=pl.len(), n_runs=pl.col("run_id").n_unique(), n_scenarios=pl.struct(keys).n_unique()
    ).collect()
    if not ((blocks["n"] == RUNS_PER_BLOCK) & (blocks["n_runs"] == RUNS_PER_BLOCK)
            & (blocks["n_scenarios"] == 1)).all():
        raise AssertionError(f"{school.name}: blocks are not {RUNS_PER_BLOCK} unique runs of one scenario")

    sizes = (
        lf.group_by(keys)
        .agg(
            n_blocks=pl.col("block").n_unique(),
            n=(pl.col("version") == "canonical").sum(),
            n_unique_runs=pl.col("run_id").filter(pl.col("version") == "canonical").n_unique(),
        )
        .collect()
    )
    if not sizes["n_blocks"].is_in([1, 2]).all():
        raise AssertionError(f"{school.name}: scenarios with more than two blocks")
    bad = sizes.filter((pl.col("n") != RUNS_PER_BLOCK) | (pl.col("n_unique_runs") != RUNS_PER_BLOCK))
    if bad.height:
        raise AssertionError(
            f"{school.name}: {bad.height}/{sizes.height} scenarios lack exactly one canonical block "
            f"of {RUNS_PER_BLOCK} unique runs (missing input column, or tied blocks)"
        )
    return sizes


def anomalies(school: School) -> dict[str, int]:
    """Counts of the known (tolerated) data anomalies, for reporting alongside the targets."""
    lf = load(school)
    n_people = school.n_students + school.n_teachers_staff
    excess = pl.col("n_infected") - pl.col("n_infected_students") - pl.col("n_infected_teachers_staff")
    row = lf.select(
        rows=pl.len(),
        inflated_rows=(pl.col("version") == "inflated").sum(),
        total_exceeds_groups=(excess > 0).sum(),
        max_total_excess=excess.max(),
        cumulative_exceeds_population=(
            (pl.col("n_infected") > n_people)
            | (pl.col("n_infected_students") > school.n_students)
            | (pl.col("n_infected_teachers_staff") > school.n_teachers_staff)
        ).sum(),
        canonical_exceeds_population=(
            (pl.col("version") == "canonical")
            & ((pl.col("n_infected_students") > school.n_students)
               | (pl.col("n_infected_teachers_staff") > school.n_teachers_staff))
        ).sum(),
        negative_endogenous=(pl.col("n_endo") < 0).sum(),
    ).collect()
    return {k: int(v) for k, v in row.row(0, named=True).items()}


def baseline(lf: pl.LazyFrame, school: School) -> pl.LazyFrame:
    """Pin the dimensions the paper holds fixed in Figs. 3-5."""
    return lf.filter(
        pl.col("version") == "canonical",
        (pl.col("susceptibility") - school.susceptibility).abs() < 1e-6,
        pl.col("classroom_isolation") == 0,
        pl.col("contact_isolation") == 0,
        pl.col("teacher_vaccine") == 0,
        pl.col("student_vaccine") == 0,
    )


# --------------------------------------------------------------------------------------------
# 3. Summaries
# --------------------------------------------------------------------------------------------

def summarize(lf: pl.LazyFrame, by: list[str], outcome: str) -> pl.DataFrame:
    """Per-group n, median, 95th percentile (linear interpolation, as numpy) and % runs > 5%."""
    col = pl.col(OUTCOMES[outcome])
    return (
        lf.group_by(by)
        .agg(
            n=pl.len(),
            median=col.median(),
            p95=col.quantile(0.95, interpolation="linear"),
            pct_gt5=(col > SIZABLE_OUTBREAK).mean() * 100,
        )
        .with_columns(outcome=pl.lit(outcome))
        .collect()
    )


# Each figure = fixed filters + grouping keys + outcomes. Mirrors the PNAS captions.
FIGURES = {
    # Fig. 3A-D: all students onsite 5 d/wk, no testing, by introduction rate (R0 1.5 and 2.25)
    "fig3_jitter": dict(
        filters=dict(cohorting="all-onsite", testing="none"), R0=[1.5, 2.25],
        by=["R0", "intro"], outcomes=["school", "school_seedfix"],
    ),
    # Fig. 3E-F: same arm, heatmap over every R0 x introduction rate, students and teachers separately
    "fig3_heatmap": dict(
        filters=dict(cohorting="all-onsite", testing="none"), R0=None,
        by=["R0", "intro"], outcomes=["students", "teachers_staff"],
    ),
    # Fig. 5: R0 1.5, ~weekly introductions, no testing, by cohorting strategy
    "fig5_jitter": dict(
        filters=dict(intro="weekly", testing="none"), R0=[1.5],
        by=["R0", "cohorting"], outcomes=["school", "school_seedfix"],
    ),
    # Fig. 4: full cohorting x testing x R0 x introduction grid, students and teachers separately
    "fig4_heatmap": dict(
        filters={}, R0=None,
        by=["R0", "intro", "cohorting", "testing"], outcomes=["students", "teachers_staff"],
    ),
}


def figure_frame(school: School, figure: str) -> pl.LazyFrame:
    """Baseline-filtered rows for one figure (shared by targets() and the plotting code)."""
    spec = FIGURES[figure]
    lf = baseline(load(school), school)
    if spec["filters"]:
        lf = lf.filter(**spec["filters"])
    if spec["R0"] is not None:
        lf = lf.filter(pl.col("R0").is_in(spec["R0"]))
    return lf


def targets(schools: list[School]) -> pl.DataFrame:
    """Summaries for every figure x school x outcome, in one tidy table."""
    frames = []
    for school in schools:
        for figure, spec in FIGURES.items():
            lf = figure_frame(school, figure)
            for outcome in spec["outcomes"]:
                summary = summarize(lf, spec["by"], outcome)
                if summary.is_empty():
                    raise ValueError(f"{school.name}/{figure}: no rows match the figure's filters")
                frames.append(summary.with_columns(figure=pl.lit(figure), school=pl.lit(school.name)))
    keys = ["figure", "school", "outcome", "R0", "intro", "cohorting", "testing"]
    df = pl.concat(frames, how="diagonal_relaxed")
    return df.select(*keys, pl.exclude(keys)).sort(keys, nulls_last=True)


# --------------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--school", choices=[*SCHOOLS, "both"], default="both")
    p.add_argument("--overwrite", action="store_true", help="re-extract parquet from raw .csv.gz")
    p.add_argument("--no-figures", action="store_true")
    args = p.parse_args(argv)

    schools = list(SCHOOLS.values()) if args.school == "both" else [SCHOOLS[args.school]]
    for s in schools:
        print(f"[{s.name}] extract -> {extract(s, overwrite=args.overwrite)}")
        sizes = validate(s)
        print(f"[{s.name}] validated: {len(sizes)} scenarios, "
              f"runs/scenario min={sizes['n'].min()} max={sizes['n'].max()}")
        print(f"[{s.name}] known anomalies: {anomalies(s)}")

    RESULTS_DIR.mkdir(exist_ok=True)
    out = RESULTS_DIR / "zenodo_targets.csv"
    targets(schools).write_csv(out, float_precision=4)
    print(f"targets -> {out}")

    if not args.no_figures:
        from . import plots

        for s in schools:
            for path in plots.make_all(s):
                print(f"figure -> {path}")


if __name__ == "__main__":
    main()
