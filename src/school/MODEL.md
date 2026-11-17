# School model: equations and parameters

The working configuration of the primary-school SEIRS+ model (Part 1B). It is what `Scenario()` and
`params.PRIMARY` give by default. Values come from the PNAS SI Appendix (`data/paper_appendix.pdf`)
unless marked **data-chosen** (picked by matching the authors' Zenodo runs) or **calibrated**.

Validated for the primary school only. The secondary school, testing and vaccination are not done yet.

| Code | Contents |
|---|---|
| `network.py` | contact network and its calibration |
| `params.py` | parameters and per-person draws |
| `run.py` | daily simulation loop and the data-chosen options |

---

## 1. Infection rate

The model is a continuous-time (Gillespie) network model. The rate at which a susceptible person
*i* becomes exposed is

$$
\text{rate}_i \;=\; \alpha_i \left[\; p\,\bar\beta\,\frac{I}{N} \;+\; \frac{1-p}{D_i^{\text{onsite}}}
\sum_{j\,\in\,\text{today's contacts of } i} \beta_j\,\delta_{ij}\,\mathbb{1}[\,j \text{ infectious}\,] \;\right]
$$

| Symbol | Meaning |
|---|---|
| $\alpha_i$ | susceptibility of *i* |
| $p$ | share of transmission from school-wide (well-mixed) contact |
| $\bar\beta$ | mean transmissibility over everyone |
| $I/N$ | fraction of the school currently infectious (presymptomatic, symptomatic or asymptomatic) |
| $D_i^{\text{onsite}}$ | *i*'s number of contacts on the everyone-onsite network, held fixed every day (**data-chosen**) |
| $\beta_j$ | transmissibility of the infectious contact *j* (`beta_pairwise_mode="infected"`) |
| $\delta_{ij}$ | degree weighting, defined below |

Source: `seirsplus_k12/models.py:2317`.

### Per-person transmissibility

$$
\beta_i \;=\; \frac{R_{0,i}}{T^{\text{pre}}_i + T^{\text{sym}}_i}
$$

$R_{0,i}$ is the individual reproduction number. The denominator is the person's total infectious period
(presymptomatic plus symptomatic or asymptomatic).

### Degree weighting

$$
\delta_i \;=\; \max\!\left(0,\; \frac{\log D_i}{\log \overline{D}}\right), \qquad
\delta_{ij} \;=\; \tfrac12\left(\delta_i + \delta_j\right)
$$

$D_i$ and $\overline{D}$ are computed from **each day's** network (**data-chosen**, see section 4).
People with no contacts get $\delta_i = 0$.

---

## 2. Parameters

### Transmission

| Parameter | Value | Source |
|---|---|---|
| $R_0$ (mean of $R_{0,i}$) | 1.5 (primary baseline) | SI |
| $R_{0,i}$ | gamma, mean $R_0$, CV 2.0 ("80/20 rule") | SI; `tests/test_params.py` matches its tables |
| $p$ (school-wide mixing) | 0.2 | SI |
| $\alpha$ (susceptibility) | 0.6 primary students, 1.0 everyone else | SI |
| Presymptomatic / asymptomatic transmission | 100% of symptomatic | SI |

### Disease progression

All residence times are gamma-distributed, drawn independently per person and per state.

| State | Mean (days) | CV | Source |
|---|---|---|---|
| Latent (E) | 3.0 | 0.6 | SI, Zenodo |
| Presymptomatic | 2.2 | 0.5 | SI, Zenodo |
| Symptomatic / asymptomatic | 4.0 | 0.4 | SI, Zenodo |

| Parameter | Value | Source |
|---|---|---|
| Asymptomatic share | 40% of primary students, 30% of everyone else (fixed at the start) | SI |

### Isolation

| Parameter | Value | Source |
|---|---|---|
| Symptomatic self-isolation | 20% of people, 1 day after symptom onset | SI |
| Isolation length | 10 days | SI |
| Isolated people | household contacts only; no school-wide mixing ($q = 0$) | SI |

### Schedule and introductions

| Parameter | Value | Source |
|---|---|---|
| Run length | 150 days | SI, Zenodo |
| Start offset | uniform 0–27 days into the 28-day schedule | Zenodo |
| Weekends | household network only | SI |
| Cohort schedules | Fig. A3 (all on site / alternate days / alternate weeks) | SI |
| Introductions | Poisson per day: 0, 1/30, 1/7 or 1 | SI |
| Who is introduced | chosen in proportion to $\alpha$ | Zenodo (students get 85–87%) |
| Initial seed | 1 exposed person, chosen uniformly | SI |

### Network

| Parameter | Value | Source |
|---|---|---|
| Students | 6 grades × 4 classes × 20, in 2 cohort blocks | SI |
| Adults | 24 teachers (one per class) + 24 staff | SI |
| `teacher_staff_degree` | 6 | **calibrated**: all 7 Zenodo network statistics within 0.7 SE |

---

## 3. Simulation loop (each day)

1. Switch to the day's network (weekend or cohort schedule) and recompute $\delta$.
2. Draw Poisson community introductions and move them S → E, restarting their latent timer.
3. Isolate yesterday's compliant symptomatic people.
4. Advance the Gillespie dynamics to the next day (step ≤ 0.1 day).

---

## 4. Choices made from the data

These differ from, or aren't settled by, the SI text. Each is an option on `run.Scenario`.

| Option | Chosen | Alternative | Evidence |
|---|---|---|---|
| `delta_network` | `"active"`: recompute $\delta$ daily (the authors' code) | `"onsite"`: fixed, as the SI text says | Cohorting, % > 5%: 10.8 / 10.3 vs Zenodo 10.0 / 9.3. Fixed gives 4.5 / 4.3 |
| `degree_denominator` | `"onsite"`: fixed contact count | `"active"`: day's contact count | Weekly intros, % > 5%: 54.6 Zenodo; `"active"` gives 65 |
| `reset_intro_timer` | `True`: introduced cases restart the latent period | `False`: the authors' code skips it | Small effect; `True` fits better |

Side effect of daily $\delta$: on household-only days $\overline{D} < 1$, so $\delta = 0$ and siblings do not
infect each other at home. School-wide mixing still runs on those days.

---

## 5. How good the fit is

Fresh validation: R0 1.5 and 2.25 × 4 introduction rates × 3 schedules (24 scenarios), new seeds,
1,000 runs each, compared with Zenodo.

| Check | Result |
|---|---|
| % of runs > 5% endogenous | within 2 SE in 19/24 scenarios, within 3 SE in 24/24 |
| Distribution shape (KS test) | p > 0.05 in 15/24 |
| Mean outbreak size | +4.7% above Zenodo (median; range +1.7% to +22%) in all 18 scenarios with introductions. **Cause not yet found** |

Fig. 5A (primary, R0 1.5, weekly introductions, % of runs > 5%):

| Schedule | Ours | Zenodo |
|---|---|---|
| All on site | 58.6 | 54.6 |
| Alternate days | 10.8 | 10.0 |
| Alternate weeks | 10.3 | 9.3 |

---

## 6. Carrying this to the vineyard (Part 2)

- **Reusable:** the structure of the rate equation, gamma residence times and overdispersed $R_{0,i}$.
- **Replace:** every number above describes COVID in schools.
- **Justify, don't copy:** $\delta$ is a hand-made weighting, and the authors' code and paper disagree on it.
- **What the Zenodo match shows:** the simulator is implemented correctly. It says nothing about whether
  these values suit vines.
