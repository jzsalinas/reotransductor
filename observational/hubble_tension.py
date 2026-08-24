"""
Hubble Tension Resolution and Dissipative Emergent Time Analysis.
Models the apparent expansion rate discrepancy between cosmic microwave background (Planck)
and local astrophysical distance ladders (SH0ES / Pantheon+ Supernovae) via the Reotransductor emergent proper time mechanism:
dtau / dt = 1 + kappa_0 * sigma(x, t).
"""

import numpy as np
from typing import Dict, Any, Tuple, Optional, List


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
        time_elapsed = max(1.0, time_elapsed)
        # Relative proper time dilation between cluster and void
        delta_tau = max(0.0, tau_cluster_mean - tau_void_mean)
        dilation_rate = delta_tau / time_elapsed

        # Predicted local Hubble constant (physically bounded between 60 and 85 km/s/Mpc)
        h0_predicted = float(np.clip(self.h0_bg * (1.0 + dilation_rate), 60.0, 85.0))
        delta_h0_pred = max(0.0, h0_predicted - self.h0_bg)

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

    def compute_3d_environmental_h0_field(
        self,
        rho_3d: np.ndarray,
        tau_3d: np.ndarray,
        tau_start_3d: Optional[np.ndarray] = None,
        scale_factor: float = 4.5,
        h0_engine: float = 0.05
    ) -> Dict[str, Any]:
        """
        Computes continuous 3D spatial field H_0(x, y, z), environmental density bins,
        and the environmental gradient dH_0 / d log10(rho / <rho>).
        """
        rho_bar = float(np.mean(rho_3d))
        rho_safe = np.maximum(1e-4, rho_3d)
        delta_rho = (rho_safe - rho_bar) / rho_bar

        tau_start = tau_start_3d if tau_start_3d is not None else np.zeros_like(tau_3d)
        delta_tau = np.maximum(0.0, tau_3d - tau_start)
        
        if scale_factor <= 1.05:
            # Recombination / homogeneous early universe: no late-time environmental gradient
            h0_field_3d = np.full_like(delta_tau, self.h0_bg)
            slope = 0.0
        else:
            time_elapsed = max(1.0, float(scale_factor - 1.0) / max(1e-4, h0_engine))
            tau_void_floor = float(np.percentile(delta_tau, 5))
            dilation_field = np.maximum(0.0, delta_tau - tau_void_floor) / time_elapsed
            h0_field_3d = np.clip(self.h0_bg * (1.0 + dilation_field), 60.0, 85.0)

            # Environmental gradient analysis: H_0 vs log10(rho / rho_bar)
            log_overdensity = np.log10(np.maximum(1e-2, rho_safe / rho_bar)).ravel()
            h0_flat = h0_field_3d.ravel()
            cov_matrix = np.cov(log_overdensity, h0_flat)
            slope = float(cov_matrix[0, 1] / max(1e-6, cov_matrix[0, 0])) if cov_matrix.shape == (2, 2) else 0.0

        # Environmental binning
        void_mask = delta_rho < -0.3
        filament_mask = (delta_rho >= -0.3) & (delta_rho < 1.0)
        group_mask = (delta_rho >= 1.0) & (delta_rho < 3.0)
        cluster_mask = delta_rho >= 3.0

        h0_void = float(np.mean(h0_field_3d[void_mask])) if np.any(void_mask) else float(np.min(h0_field_3d))
        h0_filament = float(np.mean(h0_field_3d[filament_mask])) if np.any(filament_mask) else float(np.median(h0_field_3d))
        h0_group = float(np.mean(h0_field_3d[group_mask])) if np.any(group_mask) else float(np.percentile(h0_field_3d, 75))
        h0_cluster = float(np.mean(h0_field_3d[cluster_mask])) if np.any(cluster_mask) else float(np.max(h0_field_3d))

        density_bins = np.linspace(-1.0, 1.5, 12)
        bin_centers = 0.5 * (density_bins[:-1] + density_bins[1:])
        binned_h0 = []
        binned_h0_err = []

        log_overdensity = np.log10(np.maximum(1e-2, rho_safe / rho_bar)).ravel()
        h0_flat = h0_field_3d.ravel()
        bin_idx = np.digitize(log_overdensity, density_bins) - 1
        for b in range(len(bin_centers)):
            mask = (bin_idx == b)
            if np.any(mask):
                binned_h0.append(float(np.mean(h0_flat[mask])))
                binned_h0_err.append(float(np.std(h0_flat[mask])))
            else:
                pred = self.h0_bg + slope * max(0.0, bin_centers[b] + 0.5)
                binned_h0.append(float(pred))
                binned_h0_err.append(0.3)

        return {
            "h0_void": round(h0_void, 2),
            "h0_filament": round(h0_filament, 2),
            "h0_group": round(h0_group, 2),
            "h0_cluster": round(h0_cluster, 2),
            "environmental_gradient_slope": round(slope, 3),
            "density_bin_centers": np.round(bin_centers, 2).tolist(),
            "binned_h0": np.round(binned_h0, 2).tolist(),
            "binned_h0_err": np.round(binned_h0_err, 2).tolist(),
            "h0_field_mean": round(float(np.mean(h0_field_3d)), 2),
            "h0_field_max": round(float(np.max(h0_field_3d)), 2),
            "h0_field_min": round(float(np.min(h0_field_3d)), 2)
        }

    def evaluate_engine_state(self, engine) -> Dict[str, Any]:
        """
        Evaluates the Hubble tension prediction from an active CosmologicalEngine instance.
        """
        if float(engine.scale_factor) <= 1.05:
            # Recombination (a ~ 1.0): H_0 is homogeneous baseline
            metrics = {
                "h0_background_planck": self.h0_bg,
                "h0_observed_shoes": self.h0_local_obs,
                "h0_predicted_reotransductor": self.h0_bg,
                "delta_h0_observed": round(float(self.delta_h0_obs), 2),
                "delta_h0_predicted": 0.0,
                "dilation_fraction_pct": 0.0,
                "tension_resolution_pct": 0.0,
                "is_tension_mitigated": False
            }
        else:
            tau_raw = engine.tau - engine.tau_eon_start
            tau_field = engine.to_cpu(tau_raw) if hasattr(engine, 'to_cpu') else np.asarray(tau_raw)
            rho_field = engine.to_cpu(engine.rho) if hasattr(engine, 'to_cpu') else np.asarray(engine.rho)

            # Separate cluster regions (rho > 1.2) from void regions (rho < 0.5)
            cluster_mask = rho_field > 1.2
            void_mask = rho_field < 0.5

            tau_cluster = float(np.mean(tau_field[cluster_mask])) if np.any(cluster_mask) else float(np.mean(tau_field))
            tau_void = float(np.mean(tau_field[void_mask])) if np.any(void_mask) else float(np.min(tau_field))
            time_elapsed = max(1.0, float(engine.scale_factor - 1.0) / max(1e-5, engine.H_0))

            metrics = self.predict_local_hubble_rate(tau_cluster, tau_void, time_elapsed)

        # Append continuous 3D field diagnostics
        rho_field = engine.to_cpu(engine.rho) if hasattr(engine, 'to_cpu') else np.asarray(engine.rho)
        tau_start = engine.to_cpu(engine.tau_eon_start) if hasattr(engine, 'to_cpu') else np.asarray(getattr(engine, 'tau_eon_start', None))
        env_field = self.compute_3d_environmental_h0_field(
            rho_3d=rho_field,
            tau_3d=engine.to_cpu(engine.tau) if hasattr(engine, 'to_cpu') else np.asarray(engine.tau),
            tau_start_3d=tau_start,
            scale_factor=float(engine.scale_factor),
            h0_engine=float(engine.H_0)
        )
        metrics["environmental_field"] = env_field
        return metrics
