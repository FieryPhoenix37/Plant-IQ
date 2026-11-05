import gzip

import polars as pl
import pytest

from school import zenodo
from school.zenodo import School

N_TOTAL, N_STUD, N_TS = 10, 8, 2
TINY = School("tiny", "tiny.csv.gz", "", N_TOTAL, N_STUD, N_TS, susceptibility=0.6)


def _raw_rows():
    """Rows shaped like the Zenodo file: endo = infected - introductions - 2 (taken from students)."""
    rows = []
    for rep, (infected, intro, endo_ts) in enumerate([(1, 0, 0), (2, 0, 0), (5, 1, 1), (9, 2, 2)]):
        endo = infected - intro - 2
        rows.append({
            "R0_MEAN": 1.5, "SUSCEPTIBILITY_students": 0.6,
            "AVG_INTRODUCTIONS_PER_DAY": 0.033333335,  # float32 artifact, as in the real file
            "ONSITE_CADENCE": "all-onsite", "TESTING_PROTOCOL": "none", "TESTING_CADENCE": "never",
            "TESTING_COMPLIANCE_RATE_RANDOM_students": 0.75, "TESTING_COMPLIANCE_RATE_RANDOM_teachers": 1.0,
            "ISOLATION_COMPLIANCE_RATE_POSITIVE_GROUPMATE": 0.0, "ISOLATION_COMPLIANCE_RATE_POSITIVE_CONTACT": 0.0,
            "INIT_VACCINE_UPTAKE_teachers": 0.0, "INIT_VACCINE_UPTAKE_students": 0.0, "R0_CV": 2.0,
            "VACCINE_EFFICACY": 0.9, "VACCINE_SUSCEPTIBILITY_FACTOR": 0.0,
            "VACCINE_TRANSMISSIBILITY_FACTOR": 0.0, "STOP_INTERVENTION_VACCINE_DOSE": 1,
            "runID": f"rep{rep}", "total_num_infected": infected, "total_num_introductions": intro,
            "total_num_infected_students": infected - endo_ts, "total_num_infected_teachers-staff": endo_ts,
            "total_num_introductions_students": intro, "total_num_introductions_teachers": 0,
            "peak_num_infected": infected, "degree_mean": 20.0, "degree_CV": 0.2, "cluster_coeff": 0.9,
            "assortativity": 0.1, "degree_mean_students": 19.5, "cluster_coeff_students": 0.95,
            "degree_mean_teachers-staff": 5.75,
            "total_num_endogenous_cases": endo,
            "total_num_endogenous_cases_students": endo - endo_ts,
            "total_num_endogenous_cases_teachers-staff": endo_ts,
            "total_pct_endogenous_cases": endo / N_TOTAL,
            "total_pct_endogenous_cases_students": (endo - endo_ts) / N_STUD,
            "total_pct_endogenous_cases_teachers-staff": endo_ts / N_TS,
        })
    return rows


@pytest.fixture
def tiny(tmp_path, monkeypatch):
    monkeypatch.setattr(zenodo, "DATA_DIR", tmp_path)
    monkeypatch.setattr(zenodo, "RUNS_PER_BLOCK", 4)  # the fixture is one block of 4 runs
    rows = _raw_rows()
    header = ["", "Unnamed: 0", "Unnamed: 0.1", "Unnamed: 0.1.1", *rows[0], "degree_95CI"]
    lines = [",".join(header)]
    for i, r in enumerate(rows):
        lines.append(",".join([str(i)] * 4 + [*map(str, r.values()), '"(3.0, 28.8)"']))
    with gzip.open(tmp_path / TINY.raw_file, "wt") as f:
        f.write("\n".join(lines) + "\n")
    zenodo.extract(TINY)
    return TINY


def test_extract_and_validate(tiny):
    sizes = zenodo.validate(tiny)
    assert sizes["n"].to_list() == [4]
    df = zenodo.load(tiny).collect()
    assert df["intro"].unique().to_list() == ["monthly"]
    assert df["n_endo"].to_list() == [-1, 0, 2, 5]
    # seed-corrected outcome adds back exactly one case
    assert (df["frac_endo_seedfix"] - df["frac_endo"]).to_list() == pytest.approx([1 / N_TOTAL] * 4)


def test_summary_threshold_is_strict(tiny):
    # frac_endo = [-0.1, 0.0, 0.2, 0.5]; only 0.2 and 0.5 are > 0.05
    s = zenodo.summarize(zenodo.load(tiny), ["intro"], "school").row(0, named=True)
    assert s["n"] == 4
    assert s["pct_gt5"] == pytest.approx(50.0)
    assert s["median"] == pytest.approx(0.1)


def test_threshold_boundary():
    lf = pl.LazyFrame({"g": ["a"] * 4, "frac_endo": [0.0, 0.05, 0.0500001, 0.2]})
    s = zenodo.summarize(lf, ["g"], "school").row(0, named=True)
    assert s["pct_gt5"] == pytest.approx(50.0)  # exactly 5% is not a sizable outbreak


def test_validate_rejects_broken_endo_rule(tiny, monkeypatch):
    bad = zenodo.load(tiny).collect().with_columns(n_endo=pl.col("n_endo") + 1)
    monkeypatch.setattr(zenodo, "load", lambda s: bad.lazy())
    with pytest.raises(AssertionError, match="total: .* expected \\[2\\]"):
        zenodo.validate(tiny)


def test_validate_rejects_broken_group_rule(tiny, monkeypatch):
    bad = zenodo.load(tiny).collect().with_columns(n_endo_teachers_staff=pl.col("n_endo_teachers_staff") - 1)
    monkeypatch.setattr(zenodo, "load", lambda s: bad.lazy())
    with pytest.raises(AssertionError, match="teachers/staff: .* expected \\[0\\]"):
        zenodo.validate(tiny)


def test_anomalies_counts_negative_endogenous(tiny):
    a = zenodo.anomalies(tiny)
    assert a["rows"] == 4 and a["negative_endogenous"] == 1 and a["total_exceeds_groups"] == 0


def test_baseline_filter():
    lf = pl.LazyFrame({
        "susceptibility": [0.6, 0.6000000238, 0.8, 0.6, 0.6],  # float32 0.6 must still match
        "classroom_isolation": [0.0, 0.0, 0.0, 1.0, 0.0],
        "contact_isolation": [0.0] * 5,
        "teacher_vaccine": [0.0, 0.0, 0.0, 0.0, 1.0],
        "student_vaccine": [0.0] * 5,
        "version": ["canonical"] * 5,
    })
    assert zenodo.baseline(lf, TINY).collect().height == 2


def test_label_versions_picks_lower_block(monkeypatch):
    """Two blocks of one scenario: the lower-mean block is canonical; single blocks are canonical."""
    monkeypatch.setattr(zenodo, "RUNS_PER_BLOCK", 2)
    keys = list(zenodo.SCENARIO_COLS.values())
    scen = {k: [0.0] * 6 for k in keys}
    scen["R0"] = [1.5, 1.5, 1.5, 1.5, 2.0, 2.0]  # scenario A = blocks 0 and 1, scenario B = block 2
    lf = pl.LazyFrame({**scen, "row_idx": list(range(6)), "n_infected": [10, 12, 30, 28, 5, 6]})
    out = zenodo.label_versions(lf).collect()
    assert out["block"].to_list() == [0, 0, 1, 1, 2, 2]
    assert out["version"].to_list() == ["canonical"] * 2 + ["inflated"] * 2 + ["canonical"] * 2

    tied = zenodo.label_versions(lf.with_columns(n_infected=pl.lit(7))).collect()
    assert tied.filter(R0=1.5)["version"].to_list() == ["canonical"] * 2 + ["inflated"] * 2


def test_pairwise_sign_and_significance():
    """ln(mean_col / mean_row): negative when the column group has smaller outbreaks."""
    import numpy as np

    rng = np.random.default_rng(0)
    lf = pl.LazyFrame({
        "g": ["big"] * 200 + ["small"] * 200,
        "frac_endo": np.r_[rng.uniform(0.2, 0.4, 200), rng.uniform(0.05, 0.1, 200)],
    })
    pw = zenodo.pairwise(lf, ["g"], "school")
    cell = {(r["row_g"], r["col_g"]): r for r in pw.iter_rows(named=True)}
    assert cell[("big", "small")]["log_ratio"] < 0 < cell[("small", "big")]["log_ratio"]
    assert cell[("big", "small")]["p_value"] < 0.01
    assert cell[("big", "big")]["log_ratio"] == 0 and cell[("big", "big")]["p_value"] == 1.0
