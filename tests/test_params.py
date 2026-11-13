"""Part 1B parameters: per-person draws reproduce the statistics printed in SI Tables A.2.2 and A.2.3."""

import numpy as np
import pytest

from school.params import PRIMARY, SchoolParams, individual_params

N = 200_000  # large draw so sampling error is well below the SI's printed precision


@pytest.fixture(scope="module")
def draws():
    np.random.seed(0)
    return individual_params(PRIMARY, student_ids=list(range(N // 2)), n_nodes=N)


def test_infectious_period_matches_table_a22(draws):
    period = draws["infectious_period"]
    assert period.mean() == pytest.approx(6.2, abs=0.05)  # SI: mean 6.2, std 1.9
    assert period.std() == pytest.approx(1.9, abs=0.05)


@pytest.mark.parametrize("R0, median, pct_above_1, p80", [(1.5, 0.26, 31, 2.2), (2.25, 0.39, 38, 3.3)])
def test_individual_R0_matches_table_a23(R0, median, pct_above_1, p80):
    np.random.seed(1)
    R0_i = individual_params(SchoolParams(R0=R0), student_ids=[], n_nodes=N)["R0_i"]
    assert R0_i.mean() == pytest.approx(R0, rel=0.03)
    assert np.median(R0_i) == pytest.approx(median, abs=0.01)
    assert 100 * (R0_i > 1).mean() == pytest.approx(pct_above_1, abs=1)
    assert np.percentile(R0_i, 80) == pytest.approx(p80, abs=0.1)
    top20 = np.sort(R0_i)[-N // 5:]
    assert top20.sum() / R0_i.sum() == pytest.approx(0.8, abs=0.03)  # the "80/20 rule"


def test_group_specific_values(draws):
    students, adults = slice(0, N // 2), slice(N // 2, N)
    assert set(draws["alpha"][students]) == {0.6} and set(draws["alpha"][adults]) == {1.0}
    assert set(draws["a"][students]) == {0.4} and set(draws["a"][adults]) == {0.3}
    assert np.array_equal(draws["beta"], draws["beta_asym"])  # asymptomatic transmit equally
