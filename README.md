# DS@GT Fall 2026: Plant IQ Team
A research project with Sam Brown. Fit models of how plant disease spreads (including physics-informed and epidemiological ODE models) and use them to plan management decisions such as testing, removal and replanting under budget and labor limits.

Paper that could act as a template for the work discussed:
https://www.pnas.org/doi/full/10.1073/pnas.2108909118.

## Setup

```bash
uv venv && source .venv/bin/activate
uv pip install -e ".[dev,optimize]"
pytest
```

Zenodo data (10.5281/zenodo.5237328) goes in `data/` (gitignored). See `PROJECT_PLAN.md` for the full plan.

## Layout

```
data/          Zenodo downloads (gitignored)
src/school/    Part 1: Zenodo figure targets + vendored SEIRS+ K-12 model (seirsplus_k12/)
src/vineyard/  Parts 2-3, 5: grid/kernel, simulator, policies, economics, spatial stats
src/optimize/  Part 4: policy optimization
notebooks/     figure generation
results/       cached sim outputs (parquet)
figures/
tests/
```
