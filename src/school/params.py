"""Part 1B: school-model parameters, as stated in the PNAS SI Appendix (data/paper_appendix.pdf).

Section references are to the SI. Values marked "Zenodo" were also confirmed against the per-run
summaries in the Zenodo data (see seirsplus_k12/README.md).

Per-person draws (``individual_params``) follow SI A.2.2-A.2.3:
    - residence times ~ gamma(mean, CV), drawn independently per person and per state (Table A.2.2, Zenodo);
    - individual R0 ~ gamma(R0, CV 2.0) (Table A.2.3, Zenodo);
    - transmissibility beta_i = R0_i / total infectious period_i, where the total infectious period is
      presymptomatic + (a)symptomatic (A.2.3); pre/a/symptomatic states transmit equally (Table A.2.1);
    - asymptomatic status fixed at initialization: 40% of primary students, 30% of everyone else (A.2.2);
    - susceptibility 0.6 for primary students, 1.0 for everyone else (A.2.3).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .seirsplus_k12.utilities import gamma_dist


@dataclass(frozen=True)
class SchoolParams:
    # Transmission (A.2.3, Table A.2.1)
    R0: float = 1.5
    R0_cv: float = 2.0  # "80/20 rule" overdispersion
    p_global: float = 0.2  # 20% of transmission from casual (well-mixed) contacts
    susceptibility_students: float = 0.6  # primary; secondary uses 1.0
    # Disease progression, gamma(mean days, CV) (Table A.2.2)
    latent: tuple[float, float] = (3.0, 0.6)
    presymptomatic: tuple[float, float] = (2.2, 0.5)
    symptomatic: tuple[float, float] = (4.0, 0.4)  # Isym and Iasym alike
    pct_asymptomatic_students: float = 0.4  # primary; secondary uses 0.3
    pct_asymptomatic_adults: float = 0.3
    # Isolation (A.2.6.4)
    isolation_days: float = 10
    symptomatic_isolation_compliance: float = 0.2  # with a 1-day lag after symptom onset
    # Testing (A.2.6.2), used only in testing scenarios
    testing_compliance_students: float = 0.75
    testing_compliance_teachers_staff: float = 1.0
    test_result_lag_days: int = 1
    # Simulation (A.2.6.1; Zenodo sim_duration and init_cadence_offset)
    sim_days: float = 150
    cadence_cycle_days: int = 28  # each run starts at a uniform random offset 0-27 into the cycle


PRIMARY = SchoolParams()
SECONDARY = SchoolParams(susceptibility_students=1.0, pct_asymptomatic_students=0.3)


def individual_params(params: SchoolParams, student_ids: list[int], n_nodes: int) -> dict[str, np.ndarray]:
    """Per-person model inputs for ExtSEIRSNetworkModel, drawn from numpy's global RNG like the upstream code.

    Returns rates (1/residence time) as the model expects, plus the drawn quantities for checking.
    """
    latent = gamma_dist(*params.latent, n_nodes)
    presym = gamma_dist(*params.presymptomatic, n_nodes)
    sym = gamma_dist(*params.symptomatic, n_nodes)
    R0_i = gamma_dist(params.R0, params.R0_cv, n_nodes)

    is_student = np.zeros(n_nodes, dtype=bool)
    is_student[student_ids] = True
    pct_asym = np.where(is_student, params.pct_asymptomatic_students, params.pct_asymptomatic_adults)
    alpha = np.where(is_student, params.susceptibility_students, 1.0)

    beta = R0_i / (presym + sym)
    return {
        "sigma": 1 / latent, "lamda": 1 / presym, "gamma": 1 / sym, "gamma_asym": 1 / sym,
        "beta": beta, "beta_asym": beta, "alpha": alpha, "a": pct_asym,
        "R0_i": R0_i, "infectious_period": presym + sym,
    }
