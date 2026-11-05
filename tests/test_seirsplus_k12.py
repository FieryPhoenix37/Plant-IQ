"""Smoke tests: the vendored SEIRS+ K-12 code runs on current numpy/scipy/networkx."""

import warnings

import numpy as np
import scipy.sparse

from school.seirsplus_k12.models import ExtSEIRSNetworkModel, _adj_matrix
from school.seirsplus_k12.networks import (
    generate_highschool_contact_network,
    generate_K5_school_contact_network,
)
from school.seirsplus_k12.sim_loops import run_tti_sim


def _primary(seed=0):
    np.random.seed(seed)  # upstream code uses the global numpy RNG
    return generate_K5_school_contact_network(
        num_grades=6, num_classrooms_per_grade=4, class_sizes=20, num_student_blocks=2, num_staff=24
    )


def test_adj_matrix_keeps_matrix_semantics():
    nets, _ = _primary()
    A = _adj_matrix(nets["onsite-all"])
    assert isinstance(A, scipy.sparse.csr_matrix)
    assert A.sum(axis=1).shape == (528, 1)  # 2-D like the old networkx.adj_matrix


def test_k5_network_matches_paper_sizes():
    nets, info = _primary()
    assert nets["onsite-all"].number_of_nodes() == 528
    assert (len(info["studentIDs"]), len(info["teacherIDs"]), len(info["staffIDs"])) == (480, 24, 24)
    assert {"onsite-all", "offsite-all", "onsite-block1", "onsite-block2"} <= set(nets)
    assert nets["offsite-all"].number_of_edges() < nets["onsite-block1"].number_of_edges()


def test_highschool_network_matches_paper_sizes():
    np.random.seed(0)
    nets, info = generate_highschool_contact_network()
    assert nets["onsite-all"].number_of_nodes() == 1000
    assert (len(info["studentIDs"]), len(info["teacherIDs"]), len(info["staffIDs"])) == (800, 125, 75)


def test_tti_sim_runs_with_network_cadence(capsys):
    nets, _ = _primary()
    model = ExtSEIRSNetworkModel(G=nets["onsite-all"], p=0.2, beta=0.3, sigma=1 / 3, lamda=1 / 2.2,
                                 gamma=1 / 4, initE=1, transition_mode="time_in_state", isolation_time=10)
    week = [nets["onsite-all"]] * 5 + [nets["offsite-all"]] * 2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run_tti_sim(model, T=14, max_dt=0.5, average_introductions_per_day=1 / 7,
                    testing_cadence="weekly", parameter_cadences={"G": week})
    assert model.t >= 14 or model.total_num_infected()[model.tidx] == 0
