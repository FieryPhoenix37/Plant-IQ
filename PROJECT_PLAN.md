# Project Plan: From School COVID Network Model to Vineyard Roguing Optimization

## Context

This project adapts a stochastic network epidemic model built for schools to grapevine red blotch virus (GRBV) in a vineyard. It then uses the adapted model to find vine removal and replanting strategies that beat published rules of thumb.

### Source papers and code
- **Base model:** McGee et al. 2021, *Model-driven mitigation measures for reopening schools during the COVID-19 pandemic*, PNAS 118(39): e2108909118. https://www.pnas.org/doi/full/10.1073/pnas.2108909118
  - The SI Appendix (PDF on the article page) has the network construction (§A.2.4) and parameters (§A.2.2–A.2.3).
- **Simulation framework:** SEIRS+ — https://github.com/ryansmcgee/seirsplus (`pip install seirsplus`, v1.0)
  - Relevant pieces: `ExtSEIRSNetworkModel`, `seirsplus.sim_loops.run_tti_sim`, `seirsplus.utilities.gamma_dist`, `examples/Extended_SEIRS_Workplace_TTI_Demo.ipynb`
  - The wiki pages "ExtSEIRSNetworkModel class" and "TTI Simulation Loop" document the arguments.
- **Simulation output data:** Zenodo 10.5281/zenodo.5237328 — https://zenodo.org/records/5237328
  - Files: `column_descriptions.csv`, `simData_primarySchools.csv.gz` (~3.2 GB), `simData_secondarySchools.csv.gz` (~2.0 GB). These are gzipped tab-separated files.
- **Plant disease target:** Flasco et al. 2025, *A decade of grapevine red blotch disease epidemiology reveals zonal roguing as novel disease management*, npj Viruses 3:29. https://www.nature.com/articles/s44298-025-00111-2 (PubMed 40295834)

### Caveats that shape the plan
1. **Zenodo holds simulation output only, no code.**
2. **The PNAS results probably came from an internal SEIRS+ v2.0.** The README says code is available from the author on request. Re-running on public v1.0 will be approximate. If v1.0 can't express something, document it, and consider emailing ryansmcgee@gmail.com for v2.0.
   - **Update:** the authors' school model code is public on non-default branches and is now vendored in `src/school/seirsplus_k12/` (see its README for provenance, compatibility patches, and upstream bugs kept as-is). `k12-schools` @ `3deebcf` supplies `run_tti_sim(..., parameter_cadences=...)` for weekend/cohort network switching; `v2dev` @ `f66456b` supplies `generate_K5_school_contact_network` and `generate_highschool_contact_network`. No public branch supports vaccination, so the Zenodo vaccine scenarios did come from an unreleased v2; Figs. 3–5 need none of that.
3. **Flasco et al. is a field and statistics study, not a simulation.** Its Fig. 5 roguing rules come from spatial statistics and have never been tested in simulation, so a simulation comparison is new work.

---

## Suggested repository layout

```
project/
  data/                 # Zenodo downloads (gitignored)
  src/
    school/             # Part 1B: seirsplus school replication
    vineyard/
      grid.py           # geometry, kernel
      sim.py            # vectorized stochastic simulator
      policies.py       # roguing / replanting policies (benchmarks + parametric)
      economics.py      # cost / NPV model
      stats.py          # index of dispersion, join counts, etc.
    optimize/           # Part 4
  notebooks/            # figure generation
  results/              # cached sim outputs (parquet)
  figures/
  tests/
```

**Conventions:**
- Every simulation function takes an explicit `seed` (`numpy.random.Generator`).
- Cache results to parquet.
- Use `multiprocessing` or `joblib` for replicates.

---

## Part 1 — Reproduce baseline published figures

### 1A. Exact reproduction from the Zenodo data (do this first)
1. Download `column_descriptions.csv` and identify these columns:
   - school type
   - R0
   - introduction rate
   - cohorting strategy
   - testing strategy
   - replicate ID
   - outcome (school-transmission cases as % of the population; also separate student and teacher columns, if present)
2. Stream the large files with `polars.scan_csv(..., separator="\t")` or DuckDB. Do **not** load them whole into pandas.
3. Rebuild these figures:
   - **Fig. 3A/B:** R0 = 1.5, all students in person 5 days a week, no testing, one distribution per introduction rate.
     - Draw jitter plots with the median (black) and 95th percentile (orange).
     - Print the % of runs with >5% of the population infected in school.
   - **Fig. 5A/B:** R0 = 1.5, roughly weekly introductions, no testing. Compare three cohorting arms: none, alternating days, and alternating weeks.
   - **Fig. 4:** heatmaps of the fraction of runs >5% over the testing × cohorting × (R0, introduction rate) grid.
4. Save the summary numbers to `results/zenodo_targets.csv`. They are the ground truth for 1B.

**Done when:** these panels visually match the published figures, and the ">5%" percentages match the printed values.

**Implementation notes (verified against the data; code in `src/school/zenodo.py`, run `python -m school.zenodo`):**
- Files are **comma**-separated despite the Zenodo description, with leftover pandas index columns.
- `*_pct_*` columns are **fractions** (0–1); ">5%" is `> 0.05`.
- No school-type, testing-arm-only, or replicate column as such: school = file; testing arm = `TESTING_PROTOCOL`; replicate = `runID` (1,000 unique per scenario).
- Extra scenario dimensions must be pinned to the paper's baseline: `SUSCEPTIBILITY_students` = 0.6 (primary) / 1.0 (secondary), `ISOLATION_COMPLIANCE_RATE_POSITIVE_GROUPMATE` = 0, `INIT_VACCINE_UPTAKE_*` = 0.
- Per group, in every row: endogenous = infected − introductions − 2 (students) / − 0 (teachers+staff). The seed case isn't in `introductions`, so one case too many is subtracted. Targets use the published column; `*_seedfix` variants add the case back.
- `total_num_infected` exceeds students + teachers/staff by 0–9 infections in ~15% of runs (more with more introductions).
- Raw files are fetched and MD5-verified by `zenodo.download()`; data-quality evidence is written by `school/diagnostics.py` to `results/zenodo_anomalies.csv` and `results/zenodo_block_comparison.csv`.
- Secondary school-wide % divides by **925** (800 students + 125 teachers), while the teachers+staff % divides by 200.
- Rows are blocks of 1,000 runs (`Unnamed: 0 // 1000`). **2,546 of 2,880 secondary scenarios have two blocks** with the same inputs but **different contact networks**: per-block clustering ≈ 0.161 vs ≈ 0.129 (mean degree ≈ 24.1 vs ≈ 25.8), no overlap. The high-clustering version is present in all 2,880 scenarios and never exceeds the population, so it is "canonical" and used for targets. The "alt_network" version has larger outbreaks (median ×1.02 at R0 0.75 up to ×1.31 at R0 2.5; max ×2.85), a slightly higher share of runs > 5% (+0.8 to +2.6 points; significant in 21% of scenarios), and all 72,526 counts above the population. Evidence: `results/zenodo_block_comparison.csv`. The alt network is a lead for Part 1B (e.g. the high-school generator's inter-grade-edge bug).
- **Decision (2026-11-11): use canonical only; alt_network is excluded from all targets and calibration.**
  - **Paper match:** the printed Fig. 3B (secondary, R0 1.5, daily introductions) has a median of about 40%. Canonical gives 39.6%; alt_network gives 73.6%.
  - **Method:** the paper reports "1,000 replicates for each parameter set", so the blocks were never pooled.
  - **Physical validity:** canonical is complete (all 2,880 scenarios) and never exceeds the population.
  - **Probable origin of alt_network:** it looks like a second, incomplete run. Its blocks are interleaved through the file, and it is missing 334 scenarios, all at susceptibility 0.8 with R0 ≥ 2.0.
  - **Still open:** the network-statistics cross-check in SI Appendix section A.2.4, where clustering ≈ 0.16 would confirm canonical. Not yet read.

**Part 1A results (canonical, published outcome column), % of runs > 5%:**

| | single | monthly | weekly | daily |
|---|---|---|---|---|
| Primary, R0 1.5 | 0.7 | 6.1 | 52.8 | 100 |
| Primary, R0 2.25 | 5.2 | 35.4 | 90.9 | 100 |
| Secondary, R0 1.5 | 5.9 | 35.8 | 85.5 | 100 |
| Secondary, R0 2.25 | 24.0 | 70.6 | 99.2 | 100 |

Fig. 5 (R0 1.5, weekly): primary 52.8 / 8.4 / 8.3, secondary 85.5 / 37.6 / 26.2 (none / alternating days / alternating weeks). Still to do: compare with the values printed in the published figures.

### 1B. Re-run with the authors' SEIRS+ K-12 code (needed so there is a working simulator to modify)
Start with the **primary school, Fig. 5A** condition only. Use `src/school/seirsplus_k12` (not pip v1.0.9).

Already established:
- Disease periods (verified against Zenodo per-run summaries): latent gamma(3.0, CV 0.6), presymptomatic gamma(2.2, CV 0.5), (a)symptomatic gamma(4.0, CV 0.4); individual R0 gamma(mean `R0_MEAN`, CV 2.0).
- Network sizes from the generators match the paper (528 primary, 1,000 secondary nodes).
- **Primary network calibrated (`src/school/network.py`):**
  - `teacher_staff_degree=6` (the upstream default is 10) reproduces all 7 Zenodo network statistics. Over 200 networks, every |z| is below 0.7 and the run-to-run SDs match.
  - The Zenodo `*_students` and `*_teachers-staff` statistics are measured on the single-group layers (`students-only`, `teacherstaff-only`), with self-loops excluded.
  - The statistics are identical across cohorting scenarios, so they are measured on the `onsite-all` network.
  - Both RNGs must be seeded: numpy for the generator, and the stdlib `random` module for FARZ.
  - Still to do: the secondary generator. This also tests whether the published runs used the buggy high-school generator.
- **Facts from the raw Zenodo columns** (first 200k primary rows):
  - Simulation length is 150 days. `init_cadence_offset` is a uniform 0–27-day start offset into the 28-day testing/cohorting cycle.
  - About 43 weekend days per person, which matches 150 × 2/7.
  - No symptomatic testing, contact tracing, hospitalizations or deaths.
  - Isolations still happen without testing: 11.8 per run, about 9.8 days each. This is symptomatic self-isolation.
- **SI Appendix read (`data/paper_appendix.pdf`, gitignored). Parameters are in `src/school/params.py`, with tests against Tables A.2.2–A.2.3:**
  - Asymptomatic: 40% of primary students, 30% of everyone else. Presymptomatic and asymptomatic people transmit as much as symptomatic ones.
  - β_i = R0_i / (presymptomatic + symptomatic period).
  - Draws reproduce the SI statistics: infectious period mean 6.2 / std 1.9. Individual R0 median 0.26 (R0 1.5) and 0.39 (R0 2.25); 31% / 38% above 1. The top 20% carry 80% of transmission.
  - Introductions: Poisson per day, assigned to people in proportion to susceptibility. Introductions to non-susceptible people are no-ops.
  - Quarantine network = household edges only, with no global transmission. The weekend network is the same, but global transmission stays at 20%.
  - Symptomatic self-isolation: 20%, with a 1-day lag. Testing (Mon or Mon+Thu): 1-day result lag; E 0% sensitivity, presymptomatic 75% for the first 2 days then 80%.
  - The vendored default false-negative table gives 78% after day 2, not 80%. This doesn't matter for the no-testing scenarios.
- **SI vs our network:**
  - Table A.2.4a lists teacher/staff mean degree "5", but that is the FARZ `m` input.
  - Table A.2.4b's own std 4.8 and CV² 0.8 imply a mean of 4.8/√0.8 = 5.4.
  - m = 6 reproduces A.2.4b (std 4.87, CV² 0.81, clustering 0.49, assortativity 0.44) and Zenodo (mean 5.45, std 4.84, clustering 0.486). m = 5 gives mean 4.48, std 4.14, clustering 0.44. **Keep m = 6.**
- **Runner requirements found in code review:**
  1. δ(ji) = (log D_i + log D_j) / (2 log D̄) is fixed from the `onsite-all` network (SI A.2.4.5). `update_parameters()` recomputes δ from whichever network is active, so pass `delta` explicitly, with `delta_pairwise_mode="mean"`.
  2. Open: the local-transmission denominator |C_i| comes from the active network in the code. On weekends that is household size, which concentrates pressure on siblings. The SI doesn't say whether this was also fixed. Test both against Zenodo.
- Plan: fit any remaining unknowns only on the no-intervention single-introduction scenarios. Hold out the Fig. 5A cohorting arms as validation.
- **Runner (`src/school/run.py`):** our own daily loop around the unchanged model, following SI A.2.6.1. It differs from the public `run_tti_sim` in four ways:
  - susceptibility-weighted introductions (Zenodo students' share is 85–87%, matching weighted selection at 85.7%, not uniform at 90.9%);
  - it keeps running after the school is temporarily infection-free;
  - a random 0–27-day schedule offset;
  - it counts successful introductions.
  - The seed is uniform over people (Zenodo: 90.9% students at every susceptibility level).
  - Bug found and fixed: the "skip idle day" check read stored counts, so cases introduced into an infection-free school never progressed. A regression test is in `tests/test_run.py`.
- **Calibration** (primary, R0 1.5, all on site, no testing; 1,000 runs each; outcome = Zenodo `*_seedfix`):

  | Introductions | % > 5% sim / Zenodo | Mean % sim / Zenodo | Verdict |
  |---|---|---|---|
  | single | 1.0 / 0.9 | 0.38 / 0.33 | match (KS p 0.91) |
  | monthly | 8.7 / 6.8 | 1.97 / 1.77 | borderline (z 1.6; MWU p 0.03) |
  | weekly | 57.1 / 54.6 | 6.21 / 6.03 | match (KS p 0.29) |
  | daily | 100 / 100 | 22.4 / 21.4 | sim slightly high (z 5; about 5% relative) |

  This uses the SI-consistent variant: a fixed onsite degree denominator and timer reset for introduced cases. The code-default variant ("active" denominator, no reset) overshoots, e.g. 69% at weekly.
- **Hold-out, Fig. 5A cohorting (weekly introductions): NOT reproduced yet.**

  | % > 5% | onsite denominator | active denominator | Zenodo |
  |---|---|---|---|
  | alternate days | 4.5 | 35.7 | 10.0 |
  | alternate weeks | 4.3 | 31.3 | 9.3 |

  The two denominators bracket Zenodo, so the authors' cohort-day transmission sits between "per-contact weight fixed" and "per-person pressure fixed". Next: test intermediate mechanisms against the cohorting arms, e.g. the cohort split, teacher edges, and how off-site students mix.

1. **Network** (SI §A.2.4):
   - 480 students in classes; each class is a fully connected clique plus one teacher.
   - Sibling/household edges between students.
   - 24 teachers and 24 staff connected in a teacher/staff network.
   - Generate a new random network for each replicate.
2. **Parameters** (SI §A.2.2–A.2.3):
   - latent, presymptomatic, and infectious period distributions
   - child susceptibility 0.6 relative to adults
   - 20% of symptomatic cases self-isolate; quarantine lasts 10 days
   - individual R0 drawn from a gamma distribution (overdispersed, so 80% of transmission comes from 20% of cases), calibrated to a mean of 1.5
   - network locality p ≈ 0.2, so 80% of transmission happens along network edges
3. **Simulation:**
   - Run `run_tti_sim` for 150 days, starting on a random weekday with one initial introduction.
   - Set Poisson introductions with `average_introductions_per_day`.
   - Weekends have no close contacts outside the household.
   - Where testing applies: 75% compliance for students, 100% for staff, 1-day turnaround.
4. Run 1,000 replicates in parallel and compare against `results/zenodo_targets.csv`.
5. List everything v1.0 can't express directly and how it was approximated, such as the log-degree transmission scaling and cohorting network switching (via `checkpoints` or swapping graphs).

**Done when:** the simulated % of runs >5% falls in the same range as the published value, and the remaining discrepancies are documented.

---

## Part 2 — Rewire the network to a vineyard grid (Flasco geometry)

### Geometry
- Real vineyard: 2 ha, 44 rows of 114–195 vines (~7,729 vines), with 1.2 m spacing within rows and 2.1 m between rows.
- Simplification: a 44 × 176 rectangle with coordinates `x = pos * 1.2`, `y = row * 2.1`.

### Transmission kernel
- Don't use a plain `networkx.grid_2d_graph`. Use a distance-weighted kernel, truncated at about 6 vine lags:
  - exponential: `w(d) = exp(-d / alpha)`, or
  - power law: `w(d) = (1 + d/alpha)^(-b)`
- The anisotropic spacing then naturally makes within-row spread stronger than across-row spread.
- Calibration target: Flasco Fig. 3 (join counts) shows significant aggregation mostly at 1 lag, extending to 3–6 lags inside dense foci.

### Mapping the school model onto the vineyard

| School model | Vineyard model |
|---|---|
| Close-contact network | Local vine-to-vine kernel (vector *Spissistilus festinus* walking or jumping between vines) |
| Casual contacts (1 − p) | Long-distance dispersal within the vineyard |
| Community introductions | Background inoculum from outside the vineyard (Flasco links strain 1 to outside sources) |
| E | Latent: infected, not yet infectious |
| I_pre | Infected and infectious, asymptomatic; detectable only by PCR |
| I_sym | Symptomatic; visible in the October survey |
| I_asym | Long-term asymptomatic infection |
| R | **No recovery** (GRBV infection is permanent); reuse this state as "removed / empty position" |

### Time step
- Switch from days to an annual step, or a two-phase year:
  - vector transmission season (June–July)
  - October symptom survey
- Symptom latency is 1–3+ years.
- Horizon: 20–25 years, roughly one vineyard lifespan.

### Implementation
Write a vectorized NumPy simulator in `src/vineyard/sim.py`. Compute force of infection by convolving the infectious grid with the kernel (`scipy.signal.fftconvolve`). There are three reasons not to stay in seirsplus:
1. Removal and replanting change node state in ways v1.0 doesn't support.
2. There is no recovery.
3. Part 4 needs tens of thousands of fast evaluations.

**Validate:** with control off and the same network, the new simulator's epidemic curves should statistically match seirsplus. Write this as a test.

### Calibration (no control)
1. Start at about 4% symptomatic incidence, with one focus at one edge of the block plus scattered background infections.
2. Fit β (local transmission), α (kernel scale), the influx rate, and the latency distribution so the model matches these targets:
   - symptomatic incidence going from 3.9% (year 0) to 36.4% (year 9), with an exponential progress rate r ≈ 0.227 per year
   - index of dispersion D > 1 on 2×2 quadrats, increasing over time (roughly 1.95 → 2.95)
   - asymptomatic PCR-positive vines at about 9% of tested vines in late years
3. Use ABC (approximate Bayesian computation) or a simple grid or Bayesian-optimization fit.
4. Keep the posterior or plausible parameter set; Parts 3–5 reuse it.

---

## Part 3 — Add infected-plant removal and replacement

### Mechanics, per simulated year
1. **October survey:** each symptomatic vine is detected with probability `p_detect`. Optionally, PCR-test a chosen set of vines; this also detects asymptomatic infections, at a cost.
2. **Removal:** remove detected vines plus whatever zone the policy specifies. Removal timing is a parameter: before the next vector season (default) or after it.
3. **Replanting:** after `replant_delay`, the position becomes a healthy (S) young vine.
   - A fraction ε of replants are already infected (infected nursery stock).
   - Replants produce no yield for about 3 years.
4. **Whole-block replant:** an optional policy action, triggered by an incidence threshold.

### Parameters (expose all in a config dataclass)

| Parameter | Default / range | Notes |
|---|---|---|
| `p_detect` | 0.85 (0.7–0.95) | Visual detection sensitivity |
| symptom latency | 1–3+ yr | Creates a hidden reservoir; the core reason zonal roguing helps |
| removal timing | before / after the vector season | One extra season of spread if removal is late |
| ε (infected replants) | 0 (0–2%) | Bokil et al. 2019 show this changes elimination thresholds |
| PCR sensitivity | 0.9 | |
| PCR cost per vine | parameter | |
| removal cost, replant cost | parameters | Use literature values (Ricketts et al. 2017 AJEV; Atallah et al. 2015 AJAE) |
| yield loss when infected | parameter | Losses of $2,200–$68,548/ha over the vineyard lifespan (Ricketts 2017) |
| discount rate | 3–5% | For NPV |

### Sensitivity analysis
1. Run one-at-a-time sweeps.
2. Run a Latin-hypercube or Sobol global sensitivity analysis (`SALib`).
3. Report which parameters change the **ranking** of policies, not just outcome levels.

---

## Part 4 — Optimization

### Objective
- **Primary:** 25-year NPV per hectare, or healthy-vine-years minus λ × vines removed.
- **Also track:** final incidence, total vines removed, PCR tests used, and the probability of crossing 30% incidence.
- **Plot:** the Pareto front of vines removed versus final incidence or NPV.

### Policy families, from simple to rich
1. **(a) Parametric zone shape:**
   - radii within row, across row, and diagonal
   - incidence thresholds that switch between radius sets
   - a whole-block replant threshold
   
   This family has 4–8 parameters and is directly comparable to Flasco.
2. **(b) Risk-score policy (top priority):**
   - Remove vine j if `score_j = Σ_k K_s(d_jk) · symptomatic_k (+ a term for recently removed neighbors) > τ`.
   - Optimize τ and the scoring-kernel scale.
   - Optionally, PCR-test vines with intermediate scores and remove only those that test positive.
   - This policy adapts to local density, which fixed shapes can't.
3. **(c) Learned policy (stretch goal):** a CNN that maps the observed symptom map to removal probabilities, trained with CMA-ES on its weights or with PPO.

### Optimizer
- For (a) and (b), use Optuna (TPE) or BoTorch (noisy expected improvement); CMA-ES (`cma` package) is an alternative. A grid search is fine for 2–3 parameters.
- **Common random numbers:** evaluate every candidate on the same seed set (≥200 replicates during search).
- **Final evaluation:** use fresh seeds, 1,000 replicates, with bootstrap confidence intervals.
- **Robust optimization:** sample model parameters (β, α, latency, influx, `p_detect`) from the Part 2 posterior in each replicate. Optimize the expected objective, then also report results per scenario.

---

## Part 5 — Compare against benchmark strategies

### Flasco et al. Fig. 5 zonal roguing rules
The rules use a **star / ray pattern**: N vines in each of the 8 directions (within row, across row, and the 4 diagonals). Implement them in `policies.py`:

| Condition | Vines removed besides the diseased vine |
|---|---|
| Incidence <20%, aggregated | 2 each side within row + 1 each side across row = 6 |
| Incidence <20%, not aggregated | 1 in each of 8 directions (the 3×3 block) = 8 |
| Incidence 21–50% | 3 in each of 8 directions = 24 |
| Incidence >51% | 6 in each of 8 directions = 48 |

Implementation notes:
- **The paper is inconsistent on one rule.** For the "<20%, not aggregated" case, the Discussion text says 6 surrounding vines, but the Fig. 5 legend says 8. Use the legend's 8 and record the discrepancy.
- **"Aggregated"** is decided by the index of dispersion on 2×2 quadrats with a χ² test (D > 1, p < 0.05), matching the paper's methods.
- **Incidence is measured from observable data:** symptomatic vines only. Choose block-wide or local-window incidence, make the choice explicit, and test both.
- The authors themselves say not to maintain vineyards above 30% incidence, citing Ricketts 2017, so pair the high-incidence rules with the 30% whole-block replant rule as a variant.

### Other baselines
1. No control.
2. Conventional roguing: remove symptomatic vines only.
3. Hesler et al. 2022 spatial roguing: the symptomatic vine plus 2 within-row neighbors on each side.
4. Ricketts et al. 2017: rogue individual vines below 30% incidence; replant the whole block above it.
5. Atallah et al. 2015 spatial strategy: PCR-test the immediate neighbors of symptomatic vines; remove those that test positive.

### Outputs
1. Jitter plots of outcome distributions per strategy, styled like PNAS Figs. 3 and 5.
2. A Pareto plot (vines removed versus final incidence or NPV): benchmarks as points, optimized families as curves.
3. Heatmaps in the style of PNAS Fig. 4: rank or performance of each strategy across influx rate × local transmission rate.
4. Paired statistical comparisons on the same seeds (Wilcoxon signed-rank or Mann–Whitney U), plus a pairwise log-ratio matrix in the style of PNAS Fig. 6.
5. A short write-up of where optimization beats Fig. 5, where it doesn't, and why. Hypothesis: the biggest gains are at low incidence, where fixed shapes over-remove. Under strong external influx, no roguing radius helps much.

---

## Execution order
1. Part 1A (Zenodo figures)
2. Part 1B (one primary-school panel in seirsplus)
3. Vineyard simulator, plus the validation test against seirsplus
4. Calibration to Flasco data (Part 2)
5. Part 3 mechanics and sensitivity sweeps
6. Part 5 benchmark policies
7. Part 4 optimization
8. Final comparison figures and write-up

## Additional references
- Chan & Jeger 1994, *An analytical model of plant virus disease dynamics with roguing and replanting*, J. Appl. Ecol. 31:413–427
- Bokil et al. 2019, *Optimal control of a vectored plant disease model for a crop with continuous replanting*, J. Biol. Dyn. doi:10.1080/17513758.2019.1622808
- Atallah et al. 2015, *A plant-level, spatial, bioeconomic model of plant disease diffusion and control: grapevine leafroll disease*, Am. J. Agric. Econ. doi:10.1093/ajae/aau032
- Ricketts et al. 2017, *Mitigating the economic impact of grapevine red blotch*, Am. J. Enol. Vitic. 68:127–135
- Hesler et al. 2022, *Spatial roguing reduces the incidence of leafroll disease…*, Am. J. Enol. Vitic. 73:227–236
- Flasco et al. 2024, *Investigating the latency period of grapevine red blotch virus…*, Am. J. Enol. Vitic. 75:0750018 (use for latency parameters)
