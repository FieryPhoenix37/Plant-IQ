"""Model and management parameters (Part 3 table). Every field is meant to be swept / sampled."""

from dataclasses import dataclass


@dataclass
class VineyardConfig:
    # Geometry (Flasco et al. 2025 block, simplified to a rectangle)
    n_rows: int = 44
    vines_per_row: int = 176
    within_row_spacing_m: float = 1.2
    between_row_spacing_m: float = 2.1

    # Transmission (calibrated in Part 2)
    beta: float = 0.0  # local transmission rate
    kernel: str = "exponential"  # "exponential" | "power"
    alpha: float = 1.0  # kernel scale (m)
    kernel_power: float = 2.0  # b for the power-law kernel
    kernel_max_lag: int = 6  # truncate at ~6 vine lags
    influx_rate: float = 0.0  # background inoculum from outside the block
    latency_years_mean: float = 2.0  # symptom latency, 1-3+ yr

    # Detection and removal
    p_detect: float = 0.85
    pcr_sensitivity: float = 0.9
    remove_before_vector_season: bool = True
    replant_delay_years: int = 0
    infected_replant_frac: float = 0.0  # epsilon
    replant_no_yield_years: int = 3

    # Economics (Ricketts et al. 2017; Atallah et al. 2015)
    pcr_cost_per_vine: float = 0.0
    removal_cost_per_vine: float = 0.0
    replant_cost_per_vine: float = 0.0
    yield_value_per_vine_year: float = 0.0
    infected_yield_loss_frac: float = 0.0
    discount_rate: float = 0.04

    horizon_years: int = 25
