"""Vineyard geometry and the distance-weighted transmission kernel.

Coordinates: x = pos * within_row_spacing, y = row * between_row_spacing, so the
anisotropic spacing makes within-row spread stronger than across-row spread.
"""

import numpy as np

from .config import VineyardConfig


def build_kernel(cfg: VineyardConfig) -> np.ndarray:
    """Return a (2L+1, 2L+1) kernel over row/position lags, L = cfg.kernel_max_lag, zero at the center."""
    lag = cfg.kernel_max_lag
    rows, pos = np.mgrid[-lag : lag + 1, -lag : lag + 1]
    d = np.hypot(pos * cfg.within_row_spacing_m, rows * cfg.between_row_spacing_m)
    if cfg.kernel == "exponential":
        w = np.exp(-d / cfg.alpha)
    elif cfg.kernel == "power":
        w = (1 + d / cfg.alpha) ** (-cfg.kernel_power)
    else:
        raise ValueError(f"unknown kernel {cfg.kernel!r}")
    w[lag, lag] = 0.0
    return w
