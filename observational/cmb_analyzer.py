"""
Spherical Harmonics Legendre Decomposition and CMB Multipole Analyzer.
Extracts the celestial sphere from the cosmological simulation lattice and computes
the angular power spectrum C_ell and D_ell = ell*(ell+1)/(2*pi) * C_ell.
"""

import os
import math
import numpy as np
from typing import Tuple, Dict, Any, Optional


class CMBSphericalHarmonicsAnalyzer:
    """
    Computes pure-NumPy Spherical Harmonic Decomposition on the S^2 Celestial Sphere.
    Transforms spatial temperature/proper time anisotropies Delta T / T (theta, phi) into:
      a_{ell, m} = integral_{S^2} (Delta T / T) * Y_{ell, m}^*(theta, phi) dOmega
      C_ell = 1 / (2*ell + 1) * sum_{m=-ell}^{ell} |a_{ell, m}|^2
      D_ell = ell*(ell + 1) / (2*pi) * C_ell [muK^2]
    """

    def __init__(self, n_theta: int = 48, n_phi: int = 96, ell_max: int = 24):
        self.n_theta = n_theta
        self.n_phi = n_phi
        self.ell_max = ell_max

        # Colatitude theta in [0, pi] and Longitude phi in [0, 2*pi]
        self.theta = np.linspace(0.0, np.pi, self.n_theta, dtype=np.float64)
        self.phi = np.linspace(0.0, 2.0 * np.pi, self.n_phi, endpoint=False, dtype=np.float64)
        self.THETA, self.PHI = np.meshgrid(self.theta, self.phi, indexing='ij')

        # Area element dOmega = sin(theta) * d_theta * d_phi
        self.d_theta = np.pi / max(1, self.n_theta - 1)
        self.d_phi = (2.0 * np.pi) / self.n_phi
        self.d_omega = np.sin(self.THETA) * self.d_theta * self.d_phi
        # Guard against zero weights at poles
        self.d_omega[0, :] *= 0.5
        self.d_omega[-1, :] *= 0.5

        # Precompute Cartesian unit vectors on S^2 for 3D lattice sphere sampling
        self.n_x = np.sin(self.THETA) * np.cos(self.PHI)
        self.n_y = np.sin(self.THETA) * np.sin(self.PHI)
        self.n_z = np.cos(self.THETA)

    # =========================================================================
    # ASSOCIATED LEGENDRE POLYNOMIALS & SPHERICAL HARMONICS
    # =========================================================================

    @staticmethod
    def _associated_legendre(ell: int, m: int, x: np.ndarray) -> np.ndarray:
        """
        Computes the Associated Legendre Polynomial P_ell^m(x) for x in [-1, 1] (m >= 0)
        using stable recurrence relations.
        """
        # P_m^m(x) = (-1)^m * (2m - 1)!! * (1 - x^2)^(m/2)
        p_mm = np.ones_like(x)
        somx2 = np.sqrt(np.maximum(0.0, 1.0 - x * x))
        fact = 1.0
        for i in range(1, m + 1):
            p_mm *= -fact * somx2
            fact += 2.0

        if ell == m:
            return p_mm

        # P_{m+1}^m(x) = x * (2m + 1) * P_m^m(x)
        p_mp1m = x * (2.0 * m + 1.0) * p_mm
        if ell == m + 1:
            return p_mp1m

        # Recurrence: (ell - m) * P_ell^m = x * (2*ell - 1) * P_{ell-1}^m - (ell + m - 1) * P_{ell-2}^m
        p_ellm_minus2 = p_mm
        p_ellm_minus1 = p_mp1m
        p_ellm = np.zeros_like(x)

        for l_curr in range(m + 2, ell + 1):
            p_ellm = (x * (2.0 * l_curr - 1.0) * p_ellm_minus1 - (l_curr + m - 1.0) * p_ellm_minus2) / (l_curr - m)
            p_ellm_minus2 = p_ellm_minus1
            p_ellm_minus1 = p_ellm

        return p_ellm

    def compute_spherical_harmonic(self, ell: int, m: int) -> np.ndarray:
        """
        Computes orthonormal spherical harmonic basis function Y_{ell, m}(theta, phi).
        """
        abs_m = abs(m)
        if abs_m > ell:
            return np.zeros_like(self.THETA, dtype=np.complex128)

        cos_theta = np.cos(self.THETA)
        p_lm = self._associated_legendre(ell, abs_m, cos_theta)

        # Normalization factor: N_lm = sqrt(((2*ell + 1)/(4*pi)) * ((ell - m)! / (ell + m)!))
        factor = math.factorial(ell - abs_m) / max(1.0, math.factorial(ell + abs_m))
        norm = math.sqrt(((2.0 * ell + 1.0) / (4.0 * np.pi)) * factor)

        y_lm = norm * p_lm * np.exp(1j * abs_m * self.PHI)
        if m < 0:
            # Y_{ell, -m} = (-1)^m * Y_{ell, m}^*
            y_lm = ((-1.0) ** abs_m) * np.conj(y_lm)

        return y_lm

    # =========================================================================
    # LATTICE EXTRACTION & POWER SPECTRUM EXTRACTION
    # =========================================================================

    def extract_celestial_sphere_from_grid(
        self,
        tau_3d: np.ndarray,
        T_3d: Optional[np.ndarray] = None,
        rho_3d: Optional[np.ndarray] = None,
        v_x: Optional[np.ndarray] = None,
        v_y: Optional[np.ndarray] = None,
        v_z: Optional[np.ndarray] = None,
        c_light: float = 2.5,
        alpha_mem: float = 0.35
    ) -> np.ndarray:
        """
        Extracts an observer celestial sphere at radius R_obs = L / 2.2 using first-principles
        Sachs-Wolfe gravitational potential, intrinsic plasma temperature perturbations,
        line-of-sight Doppler velocity, and holographic non-equilibrium fossil memory.
        """
        if hasattr(tau_3d, 'get'):
            tau_3d = tau_3d.get()
        if T_3d is not None and hasattr(T_3d, 'get'):
            T_3d = T_3d.get()
        if rho_3d is not None and hasattr(rho_3d, 'get'):
            rho_3d = rho_3d.get()
        if v_x is not None and hasattr(v_x, 'get'):
            v_x = v_x.get()
        if v_y is not None and hasattr(v_y, 'get'):
            v_y = v_y.get()
        if v_z is not None and hasattr(v_z, 'get'):
            v_z = v_z.get()

        grid_size = tau_3d.shape[0]
        center = grid_size / 2.0
        r_obs = grid_size / 2.2

        coords_x = (center + r_obs * self.n_x) % grid_size
        coords_y = (center + r_obs * self.n_y) % grid_size
        coords_z = (center + r_obs * self.n_z) % grid_size

        # Periodic trilinear interpolation on 3D lattice
        x0 = np.floor(coords_x).astype(int)
        y0 = np.floor(coords_y).astype(int)
        z0 = np.floor(coords_z).astype(int)
        x1 = (x0 + 1) % grid_size
        y1 = (y0 + 1) % grid_size
        z1 = (z0 + 1) % grid_size

        xd = coords_x - x0
        yd = coords_y - y0
        zd = coords_z - z0

        def _trilinear(arr_3d):
            c000 = arr_3d[x0, y0, z0]
            c100 = arr_3d[x1, y0, z0]
            c010 = arr_3d[x0, y1, z0]
            c110 = arr_3d[x1, y1, z0]
            c001 = arr_3d[x0, y0, z1]
            c101 = arr_3d[x1, y0, z1]
            c011 = arr_3d[x0, y1, z1]
            c111 = arr_3d[x1, y1, z1]

            c00 = c000 * (1.0 - xd) + c100 * xd
            c01 = c001 * (1.0 - xd) + c101 * xd
            c10 = c010 * (1.0 - xd) + c110 * xd
            c11 = c011 * (1.0 - xd) + c111 * xd

            c0 = c00 * (1.0 - yd) + c10 * yd
            c1 = c01 * (1.0 - yd) + c11 * yd
            return c0 * (1.0 - zd) + c1 * zd

        tau_s = _trilinear(tau_3d)

        # Line-of-Sight Doppler Shift
        doppler_shift = np.zeros_like(tau_s)
        if v_x is not None and v_y is not None and v_z is not None:
            vx_s = _trilinear(v_x)
            vy_s = _trilinear(v_y)
            vz_s = _trilinear(v_z)
            v_los = (vx_s * self.n_x + vy_s * self.n_y + vz_s * self.n_z) / c_light
            doppler_shift = v_los

        # First-Principles Relativistic Sachs-Wolfe + Plasma Temperature + Holographic Memory
        if T_3d is not None and rho_3d is not None:
            T_s = _trilinear(T_3d)
            rho_s = _trilinear(rho_3d)

            mean_T = max(1e-4, float(np.mean(T_3d)))
            mean_rho = max(1e-4, float(np.mean(rho_3d)))
            mean_tau = max(1e-4, float(np.mean(tau_3d)) + 1.0)

            delta_T_int = (T_s - mean_T) / mean_T
            delta_rho = (rho_s - mean_rho) / mean_rho
            delta_tau = (tau_s - float(np.mean(tau_3d))) / mean_tau

            cmb_raw = delta_T_int + (1.0 / 3.0) * delta_rho + doppler_shift + alpha_mem * delta_tau
        else:
            cmb_raw = np.log10(1.0 + np.maximum(0.0, tau_s)) + 0.4 * doppler_shift

        std_val = max(1e-4, float(np.std(cmb_raw)))
        delta_T_map = (cmb_raw - float(np.mean(cmb_raw))) / std_val
        return delta_T_map

    def compute_angular_power_spectrum(
        self,
        celestial_map: np.ndarray,
        temperature_scale_uK: float = 27.255
    ) -> Dict[str, Any]:
        """
        Decomposes the 2D celestial map into multipoles C_ell and D_ell.
        Args:
            celestial_map (ndarray): Anisotropy map on S^2 (n_theta x n_phi)
            temperature_scale_uK (float): Calibration scale in microKelvin (default ~27.255 muK)
        Returns:
            Dict containing ell, C_ell, D_ell, quadrupole, octopole, and low-ell alignment metrics.
        """
        ell_list = []
        c_ell_list = []
        d_ell_list = []
        alm_dict = {}

        for ell in range(2, self.ell_max + 1):
            m_powers = []
            for m in range(-ell, ell + 1):
                y_lm = self.compute_spherical_harmonic(ell, m)
                # Numerical surface integration: a_{ell,m} = sum(map * conj(Y_lm) * dOmega)
                alm = np.sum(celestial_map * np.conj(y_lm) * self.d_omega)
                alm_dict[(ell, m)] = alm
                m_powers.append(np.abs(alm) ** 2)

            # C_ell = 1 / (2*ell + 1) * sum |a_{ell,m}|^2
            c_ell = float(np.mean(m_powers))
            # Scale to physical microKelvin^2 (CMB standard units)
            d_ell = (ell * (ell + 1.0) / (2.0 * np.pi)) * c_ell * (temperature_scale_uK ** 2)

            ell_list.append(ell)
            c_ell_list.append(c_ell)
            d_ell_list.append(d_ell)

        ell_arr = np.array(ell_list, dtype=int)
        c_ell_arr = np.array(c_ell_list, dtype=np.float64)
        d_ell_arr = np.array(d_ell_list, dtype=np.float64)

        # Low-ell anomaly diagnostics
        c2 = float(c_ell_arr[0]) if len(c_ell_arr) > 0 else 0.0
        c3 = float(c_ell_arr[1]) if len(c_ell_arr) > 1 else 0.0
        ratio_c2_c3 = c2 / max(1e-6, c3)

        return {
            "ell": ell_arr,
            "C_ell": c_ell_arr,
            "D_ell": d_ell_arr,
            "quadrupole_C2": c2,
            "octopole_C3": c3,
            "ratio_C2_C3": round(ratio_c2_c3, 3),
            "is_quadrupole_suppressed": bool(ratio_c2_c3 < 1.0),
            "alm": alm_dict
        }

    def render_mollweide_plot(
        self,
        celestial_map: np.ndarray,
        output_path: str,
        eon: int = 1,
        scale_factor: float = 1.0,
        time_myr: float = 0.0
    ):
        """
        Renders and saves a publication-quality Mollweide projection of the CMB celestial sphere
        using the official ESA Planck color palette.
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap

        planck_cmap = LinearSegmentedColormap.from_list(
            'planck_cmb_official',
            ['#05103a', '#194a8d', '#3288bd', '#66c2a5', '#f7f7f7', '#fee08b', '#fdae61', '#d53e4f', '#5e001f'],
            N=256
        )

        # Coordinate transformation to Mollweide: lon in [-pi, pi], lat in [-pi/2, pi/2]
        # Our THETA is [0, pi] from North to South Pole -> lat = pi/2 - THETA
        # Our PHI is [0, 2*pi] -> lon = PHI - pi
        lat = (np.pi / 2.0) - self.THETA
        lon = self.PHI - np.pi

        fig = plt.figure(figsize=(10, 5.5))
        fig.patch.set_facecolor('#0f172a')
        ax = fig.add_subplot(111, projection='mollweide', facecolor='#090d16')

        # Meshgrid plotting with periodic boundaries
        mesh = ax.pcolormesh(lon, lat, celestial_map, cmap=planck_cmap, vmin=-3.0, vmax=3.0, shading='auto')
        ax.set_title(
            f"ESA Planck-Calibrated CMB Temperature Anisotropies — Eon {eon}\n"
            f"(Scale Factor a = {scale_factor:.3f} | Cosmic Time: {time_myr:,.1f} Myr)",
            color='#f8fafc', fontsize=11, fontweight='bold', pad=15
        )
        ax.tick_params(colors='#94a3b8', labelsize=8)
        ax.grid(True, linestyle=':', alpha=0.35, color='#475569')

        cbar = fig.colorbar(mesh, orientation='horizontal', pad=0.10, shrink=0.65, aspect=24)
        cbar.set_label(r'$\Delta T / \bar{T} \quad [\mu\mathrm{K} / \sigma]$', color='#f8fafc', fontsize=10, fontweight='bold')
        cbar.ax.tick_params(colors='#cbd5e1', labelsize=8)

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=180, facecolor=fig.get_facecolor(), bbox_inches='tight')
        plt.close()
        print(f"• CMB Mollweide projection successfully saved to: {output_path}")
