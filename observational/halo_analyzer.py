"""
3D Cosmological Halo Radial Density Profile & Cusp-Core Problem Analyzer.
Extracts spherical radial profiles rho(r), enclosed mass M(<r), circular velocity V_c(r),
and logarithmic density slopes gamma(r) = d ln(rho) / d ln(r) from 3D lattice simulations.
"""

from typing import Dict, Any, Tuple, Optional
import numpy as np


class HaloRadialProfileAnalyzer:
    """
    Analyzes virialized gravitational halos in 3D cosmological fields,
    evaluating radial density distributions, rotation curves, and the central logarithmic slope.
    """

    def __init__(self, grid_size: int = 32, box_size_mpc: float = 100.0, n_shells: int = 24, r_max_mpc: Optional[float] = None):
        self.grid_size = grid_size
        self.box_size_mpc = box_size_mpc
        self.dx = box_size_mpc / grid_size
        self.n_shells = n_shells
        if r_max_mpc is not None:
            self.r_max = r_max_mpc
        else:
            self.r_max = box_size_mpc / 2.0

        self.bin_edges = np.linspace(0.0, self.r_max, self.n_shells + 1, dtype=np.float64)
        self.r_centers = 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])

    def locate_halo_center(self, rho_3d: np.ndarray, phi_3d: Optional[np.ndarray] = None) -> Tuple[int, int, int]:
        """Locates the primary halo center using potential minimum or density peak."""
        if hasattr(rho_3d, 'get'):
            rho_3d = rho_3d.get()
        if phi_3d is not None and hasattr(phi_3d, 'get'):
            phi_3d = phi_3d.get()

        if phi_3d is not None:
            # Gravitational potential well minimum
            min_idx = np.unravel_index(np.argmin(phi_3d), phi_3d.shape)
            return int(min_idx[0]), int(min_idx[1]), int(min_idx[2])
        
        # Fallback to density peak
        max_idx = np.unravel_index(np.argmax(rho_3d), rho_3d.shape)
        return int(max_idx[0]), int(max_idx[1]), int(max_idx[2])

    def compute_radial_profile(
        self,
        field_3d: np.ndarray,
        center: Tuple[int, int, int]
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Computes spherically averaged radial profile centered at specified halo coordinates.
        Returns (r_centers, radial_mean, counts_per_shell).
        """
        if hasattr(field_3d, 'get'):
            field_3d = field_3d.get()
        cx, cy, cz = center
        idx = np.arange(self.grid_size)
        dx_grid = (idx[:, None, None] - cx + self.grid_size // 2) % self.grid_size - self.grid_size // 2
        dy_grid = (idx[None, :, None] - cy + self.grid_size // 2) % self.grid_size - self.grid_size // 2
        dz_grid = (idx[None, None, :] - cz + self.grid_size // 2) % self.grid_size - self.grid_size // 2

        r_dist = np.sqrt(dx_grid**2 + dy_grid**2 + dz_grid**2) * self.dx
        bin_indices = np.digitize(r_dist.ravel(), self.bin_edges) - 1

        field_flat = field_3d.ravel()
        shell_sums = np.zeros(self.n_shells, dtype=np.float64)
        shell_counts = np.zeros(self.n_shells, dtype=np.int64)

        for b in range(self.n_shells):
            mask = (bin_indices == b)
            if np.any(mask):
                shell_sums[b] = np.sum(field_flat[mask])
                shell_counts[b] = np.count_nonzero(mask)

        counts_safe = np.maximum(1, shell_counts)
        valid_shells = (shell_counts > 0)
        if np.any(valid_shells):
            raw_means = shell_sums[valid_shells] / shell_counts[valid_shells]
            radial_mean = np.interp(self.r_centers, self.r_centers[valid_shells], raw_means)
        else:
            radial_mean = np.zeros(self.n_shells, dtype=np.float64)
        return self.r_centers.copy(), radial_mean, shell_counts

    def compute_logarithmic_slope(self, r_arr: np.ndarray, rho_arr: np.ndarray) -> np.ndarray:
        """
        Computes logarithmic density slope:
          gamma(r) = d ln(rho) / d ln(r)
        Uses centered finite differences in log space.
        """
        valid = (r_arr > 0.0) & (rho_arr > 0.0)
        r_v = r_arr[valid]
        rho_v = rho_arr[valid]

        if len(r_v) < 3:
            return np.zeros_like(r_arr)

        ln_r = np.log(r_v)
        ln_rho = np.log(rho_v)

        gamma_v = np.gradient(ln_rho, ln_r)
        
        # Interpolate back to full r_arr
        gamma_full = np.interp(r_arr, r_v, gamma_v)
        return gamma_full

    def compute_circular_velocity(
        self,
        r_arr: np.ndarray,
        rho_arr: np.ndarray,
        g_const: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes enclosed mass M(<r) and circular velocity V_c(r) = sqrt(G * M(<r) / r).
        """
        dr = np.gradient(r_arr)
        # Shell volume elements: 4 * pi * r^2 * dr
        shell_mass = 4.0 * np.pi * (r_arr**2) * rho_arr * dr
        m_enclosed = np.cumsum(shell_mass)

        # Circular orbital velocity
        r_safe = np.maximum(1e-3, r_arr)
        v_circ = np.sqrt(g_const * m_enclosed / r_safe)
        return m_enclosed, v_circ

    def fit_nfw_vs_core(
        self,
        r_arr: np.ndarray,
        rho_arr: np.ndarray
    ) -> Dict[str, Any]:
        """
        Performs Pure-NumPy non-linear least-squares fitting of NFW (Cusp) and Burkert (Core) models
        against the simulation profile, determining goodness of fit and core radius.
        """
        fit_mask = (r_arr > 0.5) & (r_arr < self.r_max * 0.7) & (rho_arr > 0.0)
        r_fit = r_arr[fit_mask]
        rho_fit = rho_arr[fit_mask]

        if len(r_fit) < 4:
            return {
                "preferred_model": "Insufficient Data",
                "chi2_nfw": 0.0,
                "chi2_core": 0.0,
                "fitted_scale_radius_nfw_mpc": 8.0,
                "fitted_core_radius_mpc": 6.0
            }

        # 1. Fit NFW: rho(r) = rho_s / [ (r/r_s) * (1 + r/r_s)^2 ]
        best_chi2_nfw = float("inf")
        best_rs = 8.0
        r_s_trials = np.linspace(1.0, 30.0, 60)
        for rs in r_s_trials:
            x = np.maximum(1e-3, r_fit / rs)
            shape = 1.0 / (x * (1.0 + x)**2)
            # Optimal linear amplitude
            amp = np.sum(rho_fit * shape) / max(1e-6, np.sum(shape**2))
            model = amp * shape
            chi2 = float(np.sum(((rho_fit - model) / np.maximum(1e-3, rho_fit))**2))
            if chi2 < best_chi2_nfw:
                best_chi2_nfw = chi2
                best_rs = float(rs)

        # 2. Fit Burkert Core: rho(r) = rho_0 / [ (1 + r/r_0) * (1 + (r/r_0)^2) ]
        best_chi2_core = float("inf")
        best_r0 = 6.0
        r_0_trials = np.linspace(1.0, 30.0, 60)
        for r0 in r_0_trials:
            x = np.maximum(1e-3, r_fit / r0)
            shape = 1.0 / ((1.0 + x) * (1.0 + x**2))
            amp = np.sum(rho_fit * shape) / max(1e-6, np.sum(shape**2))
            model = amp * shape
            chi2 = float(np.sum(((rho_fit - model) / np.maximum(1e-3, rho_fit))**2))
            if chi2 < best_chi2_core:
                best_chi2_core = chi2
                best_r0 = float(r0)

        preferred = "Reotransductor Cored Profile" if best_chi2_core <= best_chi2_nfw else "NFW Cuspy Profile"

        return {
            "chi2_nfw": round(best_chi2_nfw, 3),
            "chi2_core": round(best_chi2_core, 3),
            "preferred_model": preferred,
            "fitted_scale_radius_nfw_mpc": round(best_rs, 2),
            "fitted_core_radius_mpc": round(best_r0, 2)
        }

    def evaluate_halo_diagnostics(
        self,
        rho_3d: np.ndarray,
        phi_3d: Optional[np.ndarray] = None,
        tau_3d: Optional[np.ndarray] = None
    ) -> Dict[str, Any]:
        """Full cosmological halo diagnostics suite using adaptive bound halo overdensity profile."""
        center = self.locate_halo_center(rho_3d, phi_3d)
        r_raw_arr, rho_raw, _ = self.compute_radial_profile(rho_3d, center)

        rho_min = float(np.min(rho_raw))
        rho_max = float(np.max(rho_raw))

        # If sitting on top of large cosmic expansion background (rho_min / rho_max > 0.25)
        if (rho_max > 0) and (rho_min / rho_max > 0.25):
            d_rho = np.gradient(rho_raw)
            min_indices = np.where(d_rho >= 0)[0]
            virial_idx = min_indices[0] if len(min_indices) > 0 and min_indices[0] > 3 else len(rho_raw) - 1
            r_arr = r_raw_arr[:virial_idx + 1]
            rho_vir = rho_raw[:virial_idx + 1]
            rho_bg = float(rho_vir[-1])
            rho_bound = np.maximum(1e-5, rho_vir - rho_bg + 1e-4)
        else:
            r_arr = r_raw_arr
            rho_bound = np.maximum(1e-5, rho_raw)
            virial_idx = len(rho_raw) - 1

        gamma_arr = self.compute_logarithmic_slope(r_arr, rho_bound)
        m_enc, v_circ = self.compute_circular_velocity(r_arr, rho_bound)
        fit_results = self.fit_nfw_vs_core(r_arr, rho_bound)

        # Inner slope estimation at innermost radial bin (r -> 0)
        inner_slope_gamma0 = float(gamma_arr[0]) if len(gamma_arr) > 0 else 0.0

        is_cored = bool(inner_slope_gamma0 > -0.65 or fit_results["preferred_model"] == "Reotransductor Cored Profile")

        # Temperature/Memory radial profile if available
        tau_radial = None
        if tau_3d is not None:
            _, tau_radial_raw, _ = self.compute_radial_profile(tau_3d, center)
            tau_radial = tau_radial_raw[:virial_idx + 1]

        return {
            "halo_center": {"x": center[0], "y": center[1], "z": center[2]},
            "r_mpc": np.round(r_arr, 3).tolist(),
            "rho_radial": np.round(rho_bound, 5).tolist(),
            "log_slope_gamma": np.round(gamma_arr, 3).tolist(),
            "inner_slope_gamma0": round(inner_slope_gamma0, 3),
            "is_cored": is_cored,
            "m_enclosed": np.round(m_enc, 4).tolist(),
            "v_circular": np.round(v_circ, 4).tolist(),
            "fit_results": fit_results,
            "tau_radial": np.round(tau_radial, 4).tolist() if tau_radial is not None else None
        }
