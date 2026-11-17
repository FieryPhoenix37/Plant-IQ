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
| $\bar\beta$ | mean transmissibility. The SI calls it the "population mean transmissibility of infectious individuals"; the code takes the mean of $\beta$ over all people (`models.py:2045`) |
| $I/N$ | fraction of the school currently infectious (presymptomatic, symptomatic or asymptomatic) |
| $D_i^{\text{onsite}}$ | *i*'s number of contacts on the everyone-onsite network, held fixed every day (**data-chosen**). The SI writes $\lvert C^{(i)}_G \rvert$, the number of *i*'s contacts in network $G$, without saying which network applies on weekend and cohort days |
| $\beta_j$ | transmissibility of the infectious contact *j* (`beta_pairwise_mode="infected"`) |
| $\delta_{ij}$ | degree weighting, defined below |

Sources: SI A.1.4 (the equation as printed) and `seirsplus_k12/models.py:2317` (as implemented).

### Per-person transmissibility

$$
\beta_i \;=\; \frac{R_{0,i}}{T^{\text{pre}}_i + T^{\text{sym}}_i}
$$

$R_{0,i}$ is the individual reproduction number. The denominator is the person's total infectious period
(presymptomatic plus symptomatic or asymptomatic). Source: SI A.2.3.

### Degree weighting

The SI (A.1.4.2) defines

$$
\delta_{ij} \;=\; \frac{\log D_i + \log D_j}{2 \log \overline{D}}
\;=\; \tfrac12\left(\delta_i + \delta_j\right), \qquad \delta_i = \frac{\log D_i}{\log \overline{D}}
$$

where $\overline{D}$ is the mean degree of the network. Three implementation details are ours or the code's,
not the SI's:

- $D_i$ and $\overline{D}$ are computed from **each day's** network (**data-chosen**, see section 4).
- People with no contacts get $\delta_i = 0$ (the code does this too: `models.py:2153`).
- We clip $\delta_i$ at 0. On household-only networks $\overline{D} < 1$, so $\log \overline{D} < 0$ and the
  ratio would flip sign.

---

## 2. Parameters

### Transmission

| Parameter | Value | Source |
|---|---|---|
| $R_0$ (mean of $R_{0,i}$) | 1.5 (primary baseline) | SI |
| $R_{0,i}$ | gamma, mean $R_0$, CV 2.0 ("80/20 rule") | SI Table A.2.3; `tests/test_params.py` matches its statistics |
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
| Symptomatic self-isolation | 20% of people, 1 day after symptom onset | SI A.2.6.4 |
| Isolation length | 10 days | SI |
| Isolated people | household contacts only; no school-wide mixing ($q = 0$) | SI A.2.4.3 |

### Schedule and introductions

| Parameter | Value | Source |
|---|---|---|
| Run length | 150 days | Zenodo (`sim_duration`); not stated in the SI |
| Start offset | uniform 0–27 days into the 28-day schedule | Zenodo |
| Weekends | household network only; school-wide mixing stays at 20% | SI A.2.4.4 |
| Cohort schedules | Fig. A3 (all on site / alternate days / alternate weeks) | SI |
| Introductions | Poisson per day: 0, 1/30, 1/7 or 1 | SI A.2.5 |
| Who is introduced | chosen in proportion to $\alpha$ | SI A.2.5; Zenodo agrees (students get 85–87%, weighted choice gives 85.7%) |
| Initial seed | 1 exposed person, chosen uniformly | Zenodo (seed is a student in 90.9% of runs, the student share of the school); not stated in the SI |

### Network

| Parameter | Value | Source |
|---|---|---|
| Students | 6 grades × 4 classes × 20 (480), in 2 cohort blocks | SI A.2.4.1 |
| Adults | 24 teachers (one per class) + 24 staff | SI A.2.4.1 |
| `teacher_staff_degree` | 6 | **calibrated**: all 7 Zenodo network statistics within 0.7 SE |

---

## 3. Simulation loop (each day)

1. Switch to the day's network (weekend or cohort schedule) and recompute $\delta$.
2. Draw Poisson community introductions and move them S → E, restarting their latent timer.
3. Isolate yesterday's compliant symptomatic people.
4. Advance the Gillespie dynamics to the next day. The SI says the step is capped at "a fraction of a day";
   we use 0.1 day.

Steps 1–3 follow the order in SI A.2.6.1. Testing (its step d–e) is not implemented yet.

---

## 4. Choices made from the data

These differ from, or aren't settled by, the SI text. Each is an option on `run.Scenario`.

| Option | Chosen | Alternative | Evidence |
|---|---|---|---|
| `delta_network` | `"active"`: recompute $\delta$ from each day's network | `"onsite"`: fixed from the everyone-onsite network, as SI A.2.4.5 states | Cohorting, % > 5%: 10.8 / 10.3 vs Zenodo 10.0 / 9.3. Fixed gives 4.5 / 4.3 |
| `degree_denominator` | `"onsite"`: fixed contact count | `"active"`: day's contact count | Weekly intros, % > 5%: 54.6 Zenodo; `"active"` gives 65 |
| `reset_intro_timer` | `True`: introduced cases restart the latent period | `False`: the model never resets the state timer of a manually infected person, so they can skip the latent period | Small effect; `True` fits better |

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

The Zenodo match shows the simulator is implemented correctly. It does not show that any of the values
above suit vines. Part 2 keeps the model's structure but replaces its network and its numbers.

### What carries over and what doesn't

| Part of the model | Carry over? |
|---|---|
| Form of the rate equation (local term + well-mixed term, scaled by susceptibility) | Yes |
| Gamma-distributed residence times per state | Yes, with vine values |
| Overdispersed individual transmissibility $R_{0,i}$ | Optional; only if the vine data supports it |
| Every number in section 2 | No: all describe COVID in schools |
| School network (classes, households, cohorts, weekends) | No: replaced by the vineyard graph below |
| Degree weighting $\delta$ | Only if justified, see below |

### A new graph: vines in rows

A vineyard is a grid of vines planted in rows. Vines are much closer to their neighbours **along a row**
than to the vines **in the next row**. The study block in Flasco et al. 2025 (Methods) is:

| Property | Value |
|---|---|
| Size | 2 ha 'Cabernet franc', Napa County, planted 2008 |
| Layout | 44 rows of 114–195 vines; 7,729 vines surveyed |
| Spacing within a row | 1.2 m |
| Spacing between rows | 2.1 m |
| Simplified layout (our choice) | 44 × 176 rectangle, vine at $x = \text{pos} \times 1.2$, $y = \text{row} \times 2.1$ |

Each vine is a node. Edges connect nearby vines, weighted by a distance kernel:

$$
w(d) = e^{-d/a} \quad\text{(exponential)} \qquad\text{or}\qquad w(d) = \left(1 + \tfrac{d}{a}\right)^{-b} \quad\text{(power law)}
$$

$d$ is the distance in metres. Because a row neighbour is 1.2 m away and a cross-row neighbour 2.1 m,
any decreasing kernel gives more weight to row neighbours, without a separate rule. A plain unweighted
grid (every neighbour equal) would lose this, so it isn't used. The kernel form, $a$, $b$ and the cut-off
are our modelling choices, to be fitted. Cutting off at about 6 vines matches the furthest lag at which
Flasco et al. 2025 found significant aggregation.

The school model's pieces map across as follows:

| School model | Vineyard model |
|---|---|
| Close contacts | Nearby vines, weighted by $w(d)$; spread by the three-cornered alfalfa hopper (*Spissistilus festinus*) |
| School-wide mixing ($p$) | Long-distance spread within the vineyard |
| Community introductions | Infection from outside the vineyard ("virus influx from background sources", Flasco et al. 2025) |
| Daily step | Yearly step: vector acquisition mostly in June–July, symptom survey in early to mid-October (Flasco et al. 2025) |
| Recovery (R) | None: infection is permanent. R is reused as "removed / empty position" |

### New parameters, aggregated from the vine literature

Each school value needs a vine replacement:

| Parameter | What the papers report | Source |
|---|---|---|
| Latency (first PCR-positive to symptoms) | Symptoms either when the virus is first detected or 4–12 months later; sentinel vines planted in 2015 first tested positive in 2018 | Flasco et al. 2024, AJEV, doi:10.5344/ajev.2024.24008 |
| Vector activity to symptoms | About a three-year lag, inferred from weather correlations | Flasco et al. 2025 (Discussion) |
| Infected without symptoms | 9% (162/1,783) of vines in the five subsections in 2023 were PCR-positive without symptoms | Flasco et al. 2025 |
| Local transmission strength, kernel $a$ (and $b$), outside influx | Not reported; must be fitted | Fitted to Flasco et al. 2025 |
| Roguing and replanting costs (Part 3) | Red blotch economics: Ricketts et al. 2017. Leafroll bioeconomic model to borrow structure from: Atallah et al. 2015 | Not yet extracted |

Fitting targets from Flasco et al. 2025 (no control), 2014–2023:

- **Incidence:** symptomatic incidence rose from 3.9% in 2014 to 36.4% in 2023. For the whole vineyard
  the exponential model fitted best, with rate of disease progress $r = 0.227$. A logistic model fitted
  best in the two subsections at and next to the disease focus.
- **Clustering:** the index of dispersion $D$ on 2×2-vine quadrats ranged from 1.95 in 2014 to 2.95 in
  2022 ($D > 1$ means aggregated).
- **Spatial reach:** in most subsections, aggregation was significant at one spatial lag (adjacent
  vines). In the two subsections with the highest incidence it reached further: up to three lags in
  subsection D, and up to six lags for GRBV strain 2 in subsection E.

### The degree weighting $\delta$ needs its own justification

$\delta$ (section 1) is not a measured quantity. The SI (A.1.4.2) explains the idea behind it. Because
local pressure is divided by the number of contacts, a person with many contacts gets less exposure from
each one. The authors did not want that drop to be fully proportional ("we do not assume that the
propensity of infection decreases linearly with degree"), so $\delta$ adds some weight back for
well-connected people. The SI gives this reasoning but does not derive the log formula or fit it to data.

The SI and the model code also handle $\delta$ differently when the network changes:

- SI A.2.4.5 says $\delta$ is computed once from the everyone-onsite network and reused on weekend and
  cohort days.
- In the model code, `ExtSEIRSNetworkModel.update_parameters` (`models.py:2151`) computes $\delta$ from
  the current network's degrees whenever no fixed $\delta$ is passed in. The authors' simulation loop
  switches networks by calling `update_parameters` (`sim_loops.py:95`). So unless their run scripts
  (not public) passed a fixed $\delta$, it was recomputed on every switch.
- Our runs reproduce the published cohorting results only with the recomputed version (section 4).
  This suggests the published numbers used it.

For vines, the question is whether a vine with more neighbours should lose less than proportional
exposure per neighbour. On a regular grid almost every vine has the same number of neighbours. The
exceptions are vines at row ends, block edges and gaps left by removed vines, so $\delta$ would mostly
change how those behave. Keep it only if vector data (for example, more vector activity at block edges)
supports that. Otherwise drop it, and let the kernel weights carry the spatial structure.
