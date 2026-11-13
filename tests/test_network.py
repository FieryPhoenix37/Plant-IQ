"""Part 1B network: the calibrated primary-school generator reproduces the Zenodo network statistics."""

import numpy as np

from school.network import build_primary_school_network, network_stats

# Zenodo primary runs: (mean, run-to-run SD) per statistic (see school.network docstring).
ZENODO_PRIMARY = {
    "deg_mean": (20.065, 0.0486), "deg_cv": (0.2228, 0.00424), "clustering": (0.9075, 0.00526),
    "assortativity": (0.1983, 0.0379), "deg_mean_students": (19.504, 0.0418),
    "clustering_students": (0.9511, 0.00389), "deg_mean_teachers_staff": (5.451, 0.343),
}


def test_primary_network_matches_zenodo_statistics():
    n = 20
    runs = [network_stats(build_primary_school_network(seed)[0]) for seed in range(n)]
    for key, (mean, sd) in ZENODO_PRIMARY.items():
        z = (np.mean([r[key] for r in runs]) - mean) / (sd / np.sqrt(n))
        assert abs(z) < 4, (key, z)  # a miscalibrated setting (e.g. teacher_staff_degree=10) gives |z| > 40


def test_primary_network_is_reproducible():
    a, _ = build_primary_school_network(3)
    b, _ = build_primary_school_network(3)
    assert sorted(a["onsite-all"].edges()) == sorted(b["onsite-all"].edges())
