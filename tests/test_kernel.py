import numpy as np

from vineyard.config import VineyardConfig
from vineyard.grid import build_kernel


def test_kernel_shape_and_center():
    cfg = VineyardConfig(kernel_max_lag=6)
    k = build_kernel(cfg)
    assert k.shape == (13, 13)
    assert k[6, 6] == 0.0


def test_within_row_stronger_than_across_row():
    k = build_kernel(VineyardConfig())
    c = 6
    assert k[c, c + 1] > k[c + 1, c]
    assert np.allclose(k, k[::-1, ::-1])
