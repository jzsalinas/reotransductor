"""
Hubble Tension Resolution and Dissipative Emergent Time Analysis.
Models the apparent expansion rate discrepancy between cosmic microwave background (Planck)
and local astrophysical distance ladders (SH0ES) via the Reotransductor emergent proper time mechanism:
dtau / dt = 1 + kappa_0 * sigma(x, t).
"""

import numpy as np
from typing import Dict, Any


class HubbleTensionAnalyzer:
    """
    Theoretical and Observational Analyzer for the Cosmological Hubble Tension:
      - Planck 2018 (CMB Background / Voids): H_0 = 67.36 +/- 0.54 km/s/Mpc
      - SH0ES 2022 (Local Cepheids / Supernovae in host galaxies): H_0 = 73.04 +/- 1.04 km/s/Mpc
      - Discrepancy: Delta H_0 = +5.68 km/s/Mpc (5.0 sigma tension)
    
    Reotransductor Hypothesis:
      Local astronomical clocks in dissipative, gravitationally collapsing galactic halos
      advance faster by Delta tau / t = kappa_0 * <sigma_local>, naturally producing an apparent
      higher local expansion rate (H_0_local > H_0_CMB) without modifying the global matter budget.
    """

    # Observational baselines
    H0_PLANCK_2018 = 67.36      # km/s/Mpc
    H0_PLANCK_ERR = 0.54
    H0_SHOES_2022 = 73.04       # km/s/Mpc
    H0_SHOES_ERR = 1.04

    def __init__(self, h0_bg: float = H0_PLANCK_2018, h0_local_obs: float = H0_SHOES_2022):
        self.h0_bg = float(h0_bg)
        self.h0_local_obs = float(h0_local_obs)
        self.delta_h0_obs = self.h0_local_obs - self.h0_bg
        self.fractional_tension = self.delta_h0_obs / self.h0_bg

    def predict_local_hubble_rate(
        self,
        tau_cluster_mean: float,
        tau_void_mean: float,
        time_elapsed: float
    ) -> Dict[str, Any]:
        """
        Computes the effective local Hubble constant predicted by proper time divergence.
        Args:
            tau_cluster_mean (float): Mean proper time inside virialized mass cores.
            tau_void_mean (float): Mean proper time inside cosmic voids.
            time_elapsed (float): Total coordinate time interval.
        Returns:
            Dict containing predicted H_0_local, Delta H_0, and tension resolution percentage.
        """
        time_elapsed = max(1e-4, time_elapsed)
        # Relative proper time dilation between cluster and void
        delta_tau = max(0.0, tau_cluster_mean - tau_void_mean)
        dilation_rate = delta_tau / time_elapsed

        # Predicted local Hubble constant
        h0_predicted = self.h0_bg * (1.0 + dilation_rate)
        delta_h0_pred = h0_predicted - self.h0_bg

        # Tension resolution accuracy: 100% means exact match to SH0ES 73.04 km/s/Mpc
        resolution_pct = min(100.0, (delta_h0_pred / max(1e-4, self.delta_h0_obs)) * 100.0)

        return {
            "h0_background_planck": self.h0_bg,
            "h0_observed_shoes": self.h0_local_obs,
            "h0_predicted_reotransductor": round(float(h0_predicted), 2),
            "delta_h0_observed": round(float(self.delta_h0_obs), 2),
            "delta_h0_predicted": round(float(delta_h0_pred), 2),
            "dilation_fraction_pct": round(float(dilation_rate * 100.0), 3),
            "tension_resolution_pct": round(float(resolution_pct), 1),
            "is_tension_mitigated": bool(resolution_pct >= 50.0)
        }

    def evaluate_engine_state(self, engine) -> Dict[str, Any]:
        """
        Evaluates the Hubble tension prediction from an active CosmologicalEngine instance.
        """
        tau_field = engine.tau - engine.tau_eon_start
        rho_field = engine.rho

        # Separate cluster regions (rho > 1.2) from void regions (rho < 0.5)
        cluster_mask = rho_field > 1.2
        void_mask = rho_field < 0.5

        tau_cluster = float(np.mean(tau_field[cluster_mask])) if np.any(cluster_mask) else float(np.mean(tau_field))
        tau_void = float(np.mean(tau_field[void_mask])) if np.any(void_mask) else float(np.min(tau_field))
        time_elapsed = float(engine.scale_factor - 1.0) / max(1e-5, engine.H_0)

        return self.predict_local_hubble_rate(tau_cluster, tau_void, time_elapsed)
