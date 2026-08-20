"""
3D Spatial Correlation Function and Baryon Acoustic Oscillations (BAO) Analyzer.
Computes the two-point spatial auto-correlation function xi(r) = <delta(x) delta(x+r)>
on a 3D periodic cosmological lattice using the Fourier Wiener-Khinchin theorem.
"""

from typing import Dict, Any, Tuple, Optional
import numpy as np


class BAOSpatialCorrelationAnalyzer:
    """
    Computes the isotropic 3D Two-Point Spatial Correlation Function xi(r)
    via 3D Fast Fourier Transforms (Wiener-Khinchin theorem) and spherical radial averaging.
    """

    def __init__(self, grid_size: int = 32, box_size_mpc: float = 100.0, n_bins: int = 32):
        self.grid_size = grid_size
        self.box_size_mpc = box_size_mpc
        self.dx = box_size_mpc / grid_size
        self.n_bins = n_bins
        self.r_max = box_size_mpc / 2.0  # Maximum radius within periodic minimal image

        # Precompute 3D displacement matrix in physical Mpc using periodic boundary conditions
        idx = np.arange(grid_size)
        dx_grid = (idx[:, None, None] + grid_size // 2) % grid_size - grid_size // 2
        dy_grid = (idx[None, :, None] + grid_size // 2) % grid_size - grid_size // 2
        dz_grid = (idx[None, None, :] + grid_size // 2) % grid_size - grid_size // 2

        self.r_3d = np.sqrt(dx_grid**2 + dy_grid**2 + dz_grid**2) * self.dx
        self.bin_edges = np.linspace(0.0, self.r_max, self.n_bins + 1, dtype=np.float64)
        self.r_centers = 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])

        # Precompute bin masks for ultra-fast spherical radial averaging
        self.bin_indices = np.digitize(self.r_3d.ravel(), self.bin_edges) - 1

    def compute_3d_autocorrelation(self, field_3d: np.ndarray) -> np.ndarray:
        """
        Computes 3D periodic autocorrelation volume xi_3d(r) via Wiener-Khinchin theorem:
          xi_3d = IFFT[ |FFT(delta)|^2 ] / N_voxels
        """
        mean_val = float(np.mean(field_3d))
        if abs(mean_val) > 1e-6:
            delta = (field_3d - mean_val) / mean_val
        else:
            delta = field_3d - mean_val

        # 3D Fourier Power Spectrum
        delta_fft = np.fft.fftn(delta.astype(np.float64))
        power_spectrum_3d = np.abs(delta_fft)**2

        # Inverse FFT to spatial autocorrelation
        xi_3d = np.real(np.fft.ifftn(power_spectrum_3d)) / (self.grid_size**3)
        return xi_3d

    def spherically_average_correlation(self, xi_3d: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Spherically averages 3D correlation volume into isotropic 1D radial bins xi(r).
        Returns (r_centers, xi_1d, count_per_bin).
        """
        xi_flat = xi_3d.ravel()
        xi_sum = np.zeros(self.n_bins, dtype=np.float64)
        counts = np.zeros(self.n_bins, dtype=np.int64)

        for b in range(self.n_bins):
            mask = (self.bin_indices == b)
            if np.any(mask):
                xi_sum[b] = np.sum(xi_flat[mask])
                counts[b] = np.count_nonzero(mask)

        # Avoid zero division
        counts_safe = np.maximum(1, counts)
        xi_1d = xi_sum / counts_safe
        return self.r_centers.copy(), xi_1d, counts

    def evaluate_cosmological_fields(
        self,
        rho_3d: np.ndarray,
        tau_3d: Optional[np.ndarray] = None,
        scale_factor: float = 1.0
    ) -> Dict[str, Any]:
        """
        Evaluates full cosmological BAO metrics:
          1. Matter density auto-correlation xi_rho(r)
          2. Proper time fossil memory correlation xi_tau(r)
          3. Acoustic sound horizon peak detection r_peak
          4. Correlation length r_0 (where xi(r_0) = 1.0)
        """
        # 1. Matter Correlation
        xi_3d_rho = self.compute_3d_autocorrelation(rho_3d)
        r_arr, xi_rho, counts = self.spherically_average_correlation(xi_3d_rho)

        # Scale physical distances by cosmic expansion factor a
        r_comoving = r_arr.copy()
        r_physical = r_arr * scale_factor

        # 2. Memory Field Correlation if present
        xi_tau = None
        if tau_3d is not None:
            xi_3d_tau = self.compute_3d_autocorrelation(tau_3d)
            _, xi_tau, _ = self.spherically_average_correlation(xi_3d_tau)

        # 3. Detect Correlation Length r_0: first crossing where xi(r) drops below 1.0
        r0_val = 5.0
        for i in range(len(r_arr) - 1):
            if xi_rho[i] >= 1.0 and xi_rho[i+1] < 1.0:
                # Linear interpolation
                frac = (1.0 - xi_rho[i]) / (xi_rho[i+1] - xi_rho[i])
                r0_val = float(r_arr[i] + frac * (r_arr[i+1] - r_arr[i]))
                break

        # 4. Acoustic Peak Detection on r^2 * xi(r)
        r2_xi = (r_arr**2) * xi_rho
        # Search in the acoustic range (r > 15 Mpc)
        acoustic_mask = (r_arr > 12.0) & (r_arr < self.r_max - 5.0)
        if np.any(acoustic_mask):
            sub_r = r_arr[acoustic_mask]
            sub_r2xi = r2_xi[acoustic_mask]
            peak_idx = int(np.argmax(sub_r2xi))
            r_peak = float(sub_r[peak_idx])
            peak_amp = float(sub_r2xi[peak_idx])
        else:
            r_peak = float(r_arr[int(np.argmax(r2_xi))])
            peak_amp = float(np.max(r2_xi))

        return {
            "r_comoving_mpc": np.round(r_comoving, 3).tolist(),
            "r_physical_mpc": np.round(r_physical, 3).tolist(),
            "xi_rho": np.round(xi_rho, 6).tolist(),
            "r2_xi_rho": np.round(r2_xi, 4).tolist(),
            "xi_tau": np.round(xi_tau, 6).tolist() if xi_tau is not None else None,
            "correlation_length_r0_mpc": round(r0_val, 2),
            "bao_peak_radius_mpc": round(r_peak, 2),
            "bao_peak_amplitude": round(peak_amp, 4),
            "variance_sigma2": float(np.var((rho_3d - np.mean(rho_3d)) / np.mean(rho_3d)))
        }
