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

To match the exact library versions the committed results were made with:

```bash
uv pip install -r requirements-lock.txt && uv pip install -e . --no-deps
```

See `PROJECT_PLAN.md` for the full plan.

## Part 1A: reproduce the PNAS figures from the Zenodo data

```bash
python -m school.zenodo
```

One command does everything:
1. Downloads the two raw files from Zenodo (10.5281/zenodo.5237328, ~5.2 GB) into `data/` if missing, and verifies their MD5 checksums.
2. Extracts the needed columns to `data/zenodo_*.parquet`, about 1 minute, once (`--overwrite` redoes it).
3. Validates the data's accounting rules (stops with an error if any fail).
4. Writes `results/zenodo_targets.csv` (per-figure summaries) plus `results/zenodo_anomalies.csv` and `results/zenodo_block_comparison.csv` (evidence for the data-handling choices).
5. Draws Figs. 2C/D and 3–7 for both schools into `figures/part1a/`.

Options: `--school primary|secondary`, `--no-figures`. Data notes and results are in `PROJECT_PLAN.md`, Part 1A.

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
