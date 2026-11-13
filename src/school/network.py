"""Part 1B: primary-school contact network (PNAS SI §A.2.4), calibrated to the Zenodo network columns.

The authors' generator (``seirsplus_k12.networks.generate_K5_school_contact_network``) builds 6 grades x
4 classrooms x 20 students (class cliques with household/sibling edges), 24 teachers (one per class) and
24 staff in a teacher/staff layer. A fresh network is generated per replicate.

Calibration: the only setting not fixed by the paper is ``teacher_staff_degree``. Scanning 5-7 and 10, the
value 6 reproduces all seven per-run Zenodo network statistics. Primary school, Zenodo mean (SD) over
864k runs vs 200 generated networks; every mean is within one standard error (|z| < 0.7) and the
run-to-run SDs match too:

    statistic                   Zenodo           generated
    degree_mean                 20.065 (0.049)   20.065 (0.053)
    degree_CV                   0.2228 (0.0042)  0.2227 (0.0045)
    cluster_coeff               0.9075 (0.0053)  0.9078 (0.0054)
    assortativity               0.1983 (0.0379)  0.1980 (0.0387)
    degree_mean_students        19.504 (0.042)   19.505 (0.045)   students-only layer
    cluster_coeff_students      0.9511 (0.0039)  0.9510 (0.0042)  students-only layer
    degree_mean_teachers-staff  5.451 (0.343)    5.439 (0.354)    teacher/staff-only layer

The upstream default (10) gives a teacher/staff-layer degree of 9.0. The match requires measuring without
the weight-0 self-loops that scipy ``block_diag`` leaves on student nodes (upstream bug 6 in
``seirsplus_k12/README.md``); with them, every student's degree is +2. The model sums adjacency weights
for degree, so those loops do not affect transmission and the graphs are passed to it unchanged.
"""

from __future__ import annotations

import random
import warnings

import networkx as nx
import numpy as np

from .seirsplus_k12.networks import generate_K5_school_contact_network

PRIMARY = {
    "num_grades": 6, "num_classrooms_per_grade": 4, "class_sizes": 20,
    "num_student_blocks": 2,  # two cohorts for the alternating-schedule scenarios
    "num_staff": 24,
    "teacher_staff_degree": 6,  # calibrated to Zenodo, see module docstring
}


def build_primary_school_network(seed: int) -> tuple[dict[str, nx.Graph], dict]:
    """Networks by schedule ('onsite-all', 'offsite-all', 'onsite-block1/2', plus the single layers) and node info.

    The upstream generator draws from numpy's global RNG and FARZ (teacher/staff layer) from the stdlib
    ``random`` module, so ``seed`` reseeds both.
    """
    np.random.seed(seed)
    random.seed(seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # scipy block_diag deprecation inside the vendored code
        return generate_K5_school_contact_network(**PRIMARY)


def _degree_stats(G: nx.Graph) -> tuple[float, float, float, float]:
    G = nx.Graph(G)
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    deg = np.array([k for _, k in G.degree()], dtype=float)
    return deg.mean(), deg.std() / deg.mean(), nx.average_clustering(G), nx.degree_assortativity_coefficient(G)


def network_stats(networks: dict[str, nx.Graph]) -> dict[str, float]:
    """The Zenodo per-run network statistics (renamed as in zenodo.NETWORK_COLS), self-loops excluded."""
    deg_mean, deg_cv, clustering, assortativity = _degree_stats(networks["onsite-all"])
    deg_students, _, clustering_students, _ = _degree_stats(networks["students-only"])
    deg_teachers_staff, *_ = _degree_stats(networks["teacherstaff-only"])
    return {
        "deg_mean": deg_mean, "deg_cv": deg_cv, "clustering": clustering, "assortativity": assortativity,
        "deg_mean_students": deg_students, "clustering_students": clustering_students,
        "deg_mean_teachers_staff": deg_teachers_staff,
    }
