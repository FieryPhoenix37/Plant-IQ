"""Part 1B: simulate the primary school with the authors' SEIRS+ model, following the SI's daily loop.

One replicate (``run_replicate``) = one fresh network + per-person draws + 150 days. Each day (SI A.2.6.1):
    (a) switch to the day's contact network (weekend / cohort schedule, Fig. A3);
    (b) Poisson(intro_rate) community introductions, people drawn in proportion to susceptibility (A.2.5);
    (c) isolate compliant symptomatic people 1 day after they were seen symptomatic (A.2.6.4);
    then advance the Gillespie dynamics (``ExtSEIRSNetworkModel.run_iteration``) to the next day.
Testing (A.2.6.2) is not implemented yet: the Fig. 5A targets have no testing.

Why our own loop rather than ``sim_loops.run_tti_sim``: the public loop (i) draws introductions uniformly
instead of by susceptibility (Zenodo: students get 85-87% of introductions, which matches weighted
selection, 288/336 = 85.7%, not uniform, 90.9%); (ii) stops for good when nobody is infected or isolated,
so later introductions never happen (Zenodo runs always last 150 days); (iii) has no random start offset
into the 28-day schedule (Zenodo ``init_cadence_offset`` is uniform 0-27); (iv) does not count
introductions, which the outcome columns need. The model and its Gillespie step are used unchanged.

Model choices that the SI states but the code would not do by default:
    - delta(ji) = (log D_i + log D_j) / (2 log mean D), fixed from the 'onsite-all' network (A.2.4.5).
      ``update_parameters`` would recompute it from each day's network, so it is passed explicitly.
    - quarantine network = household edges only ('offsite-all'), no global transmission (q = 0) (A.2.4.3).

Two behaviours the SI does not settle, exposed as options so they can be tested against Zenodo:
    - ``degree_denominator``: local pressure is divided by the susceptible person's number of contacts.
      'active' (the code) uses the day's network, so on weekends a sibling's pressure is undiluted;
      'onsite' (default) keeps the everyone-onsite count, like delta.
    - ``reset_intro_timer``: the model never resets the state timer of a person infected by an
      introduction, so mid-run introductions skip their latent period. True (default) resets it, as the SI
      describes ("they become exposed"); False reproduces the public code.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import polars as pl
from joblib import Parallel, delayed

from .network import build_primary_school_network
from .params import PRIMARY, SchoolParams, individual_params
from .seirsplus_k12.models import ExtSEIRSNetworkModel, _adj_matrix

# Fig. A3 schedules: network on each weekday (Mon-Fri) of a 2-week cycle; weekends are always offsite.
COHORT_SCHEDULES = {
    "all-onsite": [["onsite-all"] * 5, ["onsite-all"] * 5],
    "alternate-blocks-daily": [["onsite-block1", "onsite-block2"] * 2 + ["onsite-block1"],
                               ["onsite-block2", "onsite-block1"] * 2 + ["onsite-block2"]],
    "alternate-blocks-weekly": [["onsite-block1"] * 5, ["onsite-block2"] * 5],
}
INTRO_RATES = {"single": 0.0, "monthly": 1 / 30, "weekly": 1 / 7, "daily": 1.0}


@dataclass(frozen=True)
class Scenario:
    intro: str = "weekly"  # key of INTRO_RATES; every run also starts with one seed case on day 0
    cohorting: str = "all-onsite"  # key of COHORT_SCHEDULES
    params: SchoolParams = PRIMARY
    degree_denominator: str = "onsite"  # 'onsite' | 'active', see module docstring
    reset_intro_timer: bool = True
    max_dt: float = 0.1  # max Gillespie step (days), so daily actions are never skipped


def _network_for_day(cohorting: str, cycle_day: int) -> str:
    weekday, week = cycle_day % 7, (cycle_day // 7) % 2
    return "offsite-all" if weekday >= 5 else COHORT_SCHEDULES[cohorting][week][weekday]


def _log_degree_delta(G) -> np.ndarray:
    # Weighted degree as the model computes it (networkx .degree() would count weight-0 self-loops as +2).
    degree = np.asarray(_adj_matrix(G).sum(axis=0), dtype=float).ravel()
    with np.errstate(divide="ignore"):
        delta = np.log(degree) / np.log(degree.mean())
    delta[np.isneginf(delta)] = 0.0
    return delta


def run_replicate(scenario: Scenario, seed: int) -> dict:
    """One 150-day run; returns Zenodo-style outcome counts (see ``school.zenodo.COUNT_COLS``)."""
    p = scenario.params
    nets, info = build_primary_school_network(seed)  # reseeds numpy and stdlib random
    G = nets["onsite-all"]
    n = G.number_of_nodes()
    students = np.zeros(n, dtype=bool)
    students[info["studentIDs"]] = True

    draws = individual_params(p, info["studentIDs"], n)
    model_draws = {k: v for k, v in draws.items() if k not in ("R0_i", "infectious_period")}
    delta = _log_degree_delta(G)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = ExtSEIRSNetworkModel(
            G=G, G_Q=nets["offsite-all"], p=p.p_global, q=0, isolation_time=p.isolation_days,
            delta=delta, delta_Q=delta, delta_pairwise_mode="mean", beta_pairwise_mode="infected",
            initE=1, transition_mode="time_in_state", **model_draws,
        )
    onsite_degree = model.degree.copy()
    model.tmax = p.sim_days

    offset = np.random.randint(p.cadence_cycle_days)
    isolation_compliant = np.random.rand(n) < p.symptomatic_isolation_compliance
    intro_weights = draws["alpha"] / draws["alpha"].sum()
    intro_rate = INTRO_RATES[scenario.intro]
    seed_is_student = bool(students[model.X.flatten() == model.E][0])

    n_intro = np.zeros(n, dtype=int)
    symptomatic_queue: list[np.ndarray] = [np.array([], dtype=int)]  # 1-day isolation lag
    infected = (model.E, model.I_pre, model.I_sym, model.I_asym, model.H,
                model.Q_E, model.Q_pre, model.Q_sym, model.Q_asym)
    inactive = (model.S, model.R, model.F)  # nothing can change: no infection, no isolation timer running
    peak = 1
    current_network = "onsite-all"

    for day in range(int(p.sim_days)):
        # (a) network for the day
        name = _network_for_day(scenario.cohorting, (offset + day) % p.cadence_cycle_days)
        if name != current_network:
            model.parameters.update({"G": nets[name]})
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model.update_parameters()
            current_network = name
        if scenario.degree_denominator == "onsite":
            model.degree = onsite_degree

        # (b) community introductions
        k = np.random.poisson(intro_rate)
        if k:
            for node in np.random.choice(n, size=k, replace=False, p=intro_weights):
                state = model.X[node, 0]
                if state in (model.S, model.Q_S):
                    model.X[node] = model.E if state == model.S else model.Q_E
                    if scenario.reset_intro_timer:
                        model.timer_state[node] = 0.0
                    n_intro[node] += 1

        # (c) isolate yesterday's compliant symptomatic people; queue today's
        for node in symptomatic_queue.pop(0):
            model.set_isolation(node, True)
        symptomatic_queue.append(np.flatnonzero((model.X.flatten() == model.I_sym) & isolation_compliant))

        # advance the dynamics to the next day; skip straight there when nothing can happen. Read the live
        # states (model.X), not model.total_num_infected(): its stored counts only refresh inside
        # run_iteration, so they miss cases just introduced above.
        while model.t < day + 1:
            if np.isin(model.X, inactive).all():
                model.timer_state += day + 1 - model.t
                model.t = float(day + 1)
                break
            with np.errstate(divide="ignore"):  # 1/0 rates of unused transitions (hospital, reinfection) = never
                model.run_iteration(max_dt=scenario.max_dt)
            peak = max(peak, int(np.isin(model.X, infected).sum()))

    ever = ~np.isin(model.X.flatten(), (model.S, model.Q_S))
    out = {
        "seed": seed, "cadence_offset": offset, "seed_is_student": seed_is_student,
        "n_infected": int(ever.sum()),
        "n_infected_students": int(ever[students].sum()),
        "n_infected_teachers_staff": int(ever[~students].sum()),
        "n_introductions": int(n_intro.sum()),
        "n_introductions_students": int(n_intro[students].sum()),
        "n_introductions_teachers": int(n_intro[~students].sum()),
        "n_peak_infected": peak,
        "n_latent_at_end": int(np.isin(model.X, (model.E, model.Q_E)).sum()),
    }
    # Endogenous = infected - introductions - the seed case (correct count; compare with Zenodo *_seedfix).
    out["n_endo"] = out["n_infected"] - out["n_introductions"] - 1
    out["frac_endo"] = out["n_endo"] / n
    out["frac_endo_students"] = (out["n_infected_students"] - out["n_introductions_students"]
                                 - seed_is_student) / int(students.sum())
    return out


def run_scenario(scenario: Scenario, n_runs: int = 1000, first_seed: int = 0, n_jobs: int = -1) -> pl.DataFrame:
    """``n_runs`` independent replicates in parallel (seeds first_seed .. first_seed + n_runs - 1)."""
    rows = Parallel(n_jobs=n_jobs)(delayed(run_replicate)(scenario, s) for s in range(first_seed, first_seed + n_runs))
    return pl.DataFrame(rows).with_columns(intro=pl.lit(scenario.intro), cohorting=pl.lit(scenario.cohorting))

