# Vendored SEIRS+ for the K-12 school model

Copied from <https://github.com/ryansmcgee/seirsplus> (MIT, see `LICENSE`). This is the authors'
model code for McGee et al. 2021; the pip release (v1.0.9) lacks the school network generators and
scheduled network switching.

| File | Source | Why |
|---|---|---|
| `models.py`, `sim_loops.py`, `utilities.py`, `FARZ.py` | branch `k12-schools` @ `3deebcf` (2020-11-09) | `run_tti_sim(..., max_dt, parameter_cadences)`: per-day network/parameter schedules (weekends, cohorting) and a max Gillespie step so daily testing/isolation still tick |
| `networks.py` | branch `v2dev` @ `f66456b` (2021-10-28) | latest `generate_K5_school_contact_network` (fixes a `block_diag` bug on `k12-schools`, reconnects isolated nodes, returns `network_info`) and the only `generate_highschool_contact_network` |

Not taken: the `v2dev` `CompartmentNetworkModel` (generic v2 rewrite published after the paper, no
log-degree scaling or school sim loop). No public branch supports the vaccination scenarios in the
Zenodo data, so those came from an unreleased internal v2.

## Local changes (compatibility only, no behavior change)

- `networks.py`: `import FARZ` -> `from . import FARZ`.
- `networks.py`, `models.py`: `networkx.adj_matrix` (removed in networkx 3) -> `_adj_matrix`, which
  returns `scipy.sparse.csr_matrix(networkx.adjacency_matrix(G))`. The cast matters: networkx 3
  returns sparse *arrays*, where `*` is elementwise and `.sum(axis=1)` is 1-D; the model relies on
  *matrix* semantics.
- `networks.py`: `networkx.from_scipy_sparse_matrix` -> `networkx.from_scipy_sparse_array`.

Verified running on numpy 2.4, scipy 1.17, networkx 3.6 (`tests/test_seirsplus_k12.py`).

## Upstream behavior kept as-is (relevant to Part 1B)

These are left unfixed so simulations match what the authors ran. Fix them only as a documented variant.

1. `generate_highschool_contact_network`: inter-grade edges use the grade-local index `i` instead of
   `i + gradeStartIdx`, so every inter-grade edge starts from a grade-0 student.
2. `generate_highschool_contact_network`: `selectedTeachers` is never appended to, so a student can
   draw the same teacher more than once (fewer than 6 distinct teachers).
3. `generate_highschool_contact_network`: the per-grade teacher slice has `+1`, overlapping the next
   grade's first teacher.
4. `ExtSEIRSNetworkModel.run_iteration`: when the step is capped at `max_dt`, isolation timers
   advance by the uncapped `tau` rather than `max_dt`, so isolation can end early.
5. `models.py` L~899-900: quarantine `delta_Q` pairwise terms transpose in the opposite order to the
   non-quarantine `delta` terms (L~873-874). Possibly intentional; unverified.
6. `generate_K5_school_contact_network`: `scipy.sparse.block_diag` keeps each classroom block's zero
   diagonal as stored entries, which networkx turns into weight-0 self-loops (one per student; 1 in
   the high-school generator). Transmission is unaffected (the model adjacency diagonal is 0), but
   networkx degree counts each loop as +2 (student degree 22.5 with, 20.5 without vs Zenodo ~19.5).
   Exclude self-loops when comparing to the Zenodo degree statistics. Whether the authors' 2020
   scipy/networkx did the same is unknown.

The Zenodo output records per-run network statistics (`degree_mean_students`, `cluster_coeff`, ...),
so whether the published runs used this exact generator code can be tested directly.

## Parameters verified against the Zenodo per-run summaries

| Parameter | Value | Zenodo check (primary, first 54k runs) |
|---|---|---|
| latent period | gamma(mean 3.0, CV 0.6) | mean 3.000, CV 0.599 |
| presymptomatic period | gamma(mean 2.2, CV 0.5) | mean 2.200, CV 0.499 |
| (a)symptomatic period | gamma(mean 4.0, CV 0.4) | mean 4.000, CV 0.399 |
| individual R0 | gamma(mean `R0_MEAN`, CV 2.0) | CV 1.985-1.989 at every R0 |

The workplace demo notebook's other values (R0 CV 0.2, p 0.4, 14-day isolation) are **not** the
school values; the paper uses 80% local transmission and 10-day isolation.
