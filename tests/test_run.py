"""Part 1B runner: reproducibility and outcome accounting of single replicates."""

import pytest

from school.run import COHORT_SCHEDULES, Scenario, _network_for_day, run_replicate


def test_replicate_is_reproducible():
    assert run_replicate(Scenario(intro="weekly"), 7) == run_replicate(Scenario(intro="weekly"), 7)


@pytest.mark.parametrize("intro", ["single", "daily"])
def test_outcome_accounting(intro):
    r = run_replicate(Scenario(intro=intro), 3)
    assert r["n_infected"] == r["n_infected_students"] + r["n_infected_teachers_staff"]
    assert r["n_introductions"] == r["n_introductions_students"] + r["n_introductions_teachers"]
    assert r["n_endo"] == r["n_infected"] - r["n_introductions"] - 1 >= 0
    assert 1 <= r["n_peak_infected"] <= r["n_infected"] <= 528
    if intro == "single":
        assert r["n_introductions"] == 0


def test_schedules_follow_fig_a3():
    week = lambda cohorting, start: [_network_for_day(cohorting, d) for d in range(start, start + 7)]
    assert week("all-onsite", 0) == ["onsite-all"] * 5 + ["offsite-all"] * 2
    assert week("alternate-blocks-daily", 0)[:5] == ["onsite-block1", "onsite-block2"] * 2 + ["onsite-block1"]
    assert week("alternate-blocks-daily", 7)[:5] == ["onsite-block2", "onsite-block1"] * 2 + ["onsite-block2"]
    assert week("alternate-blocks-weekly", 7)[:5] == ["onsite-block2"] * 5
    assert set(COHORT_SCHEDULES) == {"all-onsite", "alternate-blocks-daily", "alternate-blocks-weekly"}


def test_introductions_into_a_quiet_school_progress():
    # Regression: cases introduced while nobody was infected used to stay latent forever, because the
    # loop read the model's stored counts (refreshed only inside run_iteration) to decide to skip a day.
    r = run_replicate(Scenario(intro="monthly"), 0)
    assert r["n_introductions"] > 0
    assert r["n_latent_at_end"] <= 2  # only cases introduced in the last few days may still be latent


def test_delta_is_nonnegative_on_every_schedule_network():
    import numpy as np

    from school.network import build_primary_school_network
    from school.run import _log_degree_delta

    nets, _ = build_primary_school_network(0)
    for name in ("onsite-all", "onsite-block1", "onsite-block2", "offsite-all"):
        delta = np.clip(_log_degree_delta(nets[name]), 0, None)
        assert np.isfinite(delta).all() and (delta >= 0).all(), name
