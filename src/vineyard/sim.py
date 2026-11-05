"""Annual-step stochastic simulator (vector season June-July, October symptom survey).

States per vine: S, E (latent), I_pre (infectious, PCR-detectable only), I_sym,
I_asym, R (removed / empty position). No recovery. Force of infection is the
infectious grid convolved with the kernel via scipy.signal.fftconvolve.
"""

import numpy as np

from .config import VineyardConfig

S, E, I_PRE, I_SYM, I_ASYM, R = range(6)


def simulate(cfg: VineyardConfig, policy, seed: np.random.Generator):
    raise NotImplementedError
