"""Part 1A: stream the Zenodo simulation output (10.5281/zenodo.5237328) and rebuild Figs. 3-5 targets.

The data files are gzipped tab-separated and several GB each; always scan lazily
(polars / DuckDB), never load whole into pandas.
"""

from pathlib import Path

ZENODO_RECORD = "https://zenodo.org/records/5237328"
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
RESULTS_DIR = Path(__file__).resolve().parents[2] / "results"

FILES = {
    "columns": "column_descriptions.csv",
    "primary": "simData_primarySchools.csv.gz",
    "secondary": "simData_secondarySchools.csv.gz",
}


def scan_sims(school: str = "primary"):
    """Return a lazy polars frame over the simulation output for one school type."""
    import polars as pl

    return pl.scan_csv(DATA_DIR / FILES[school], separator="\t")


def build_targets():
    """Compute % of runs >5% infected (and medians / 95th pct) per condition -> results/zenodo_targets.csv."""
    raise NotImplementedError
