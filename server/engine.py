"""
Headless Cosmological Computation Engine for 24/7 Server Execution
Optimized for multi-core CPUs (Dell PowerEdge R820) with automated checkpointing, snapshot archiving, and safe simulation reset.
"""

import os
import json
import time
import io
import csv
import shutil
import re
import numpy as np
from server.notifier import TelegramNotifier
from server.physics_units import CosmologicalUnits, FundamentalConstants, PlanckScales

class CosmologicalEngine:
    """
    Autonomous 3D Cosmological Physics Engine.
    Executes Navier-Stokes, Poisson Gravity, Onsager Emergent Time, and Bekenstein Quantum Bounce.
    Supports unified CPU (NumPy) and GPU (CuPy / CUDA) hardware execution.
    """

    def __init__(self, grid_size=32, box_size_mpc=None, checkpoint_dir="checkpoints", auto_resume=True, force_reset=False, initial_speed=20, seed=42, use_gpu=False):
        self.grid_size = grid_size
        env_box = float(os.getenv("REOTRANSDUCTOR_BOX_SIZE_MPC", "100.0"))
        self.box_size_mpc = float(box_size_mpc) if box_size_mpc is not None else env_box
        self.initial_speed = int(initial_speed)
        self.checkpoint_dir = checkpoint_dir
        self.snapshots_dir = os.path.join(self.checkpoint_dir, "snapshots")
        os.makedirs(self.checkpoint_dir, exist_ok=True)
        os.makedirs(self.snapshots_dir, exist_ok=True)
        self.history_file = os.path.join(self.checkpoint_dir, "history.json")

        # Hardware Backend Selection (Single Source of Truth: CPU / GPU)
        self.use_gpu = bool(use_gpu)
        if self.use_gpu:
            try:
                import cupy as cp
                self.xp = cp
                self.gpu_available = True
            except Exception:
                self.xp = np
                self.use_gpu = False
                self.gpu_available = False
        else:
            self.xp = np
            self.gpu_available = False

        # Centralized Seed & Reproducible Pseudo-Random Number Generator
        self.seed = int(seed)
        self.rng = np.random.default_rng(self.seed)

        # Telegram Notifier (loads from local gitignored telegram_config.json)
        self.notifier = TelegramNotifier()

        # Cosmological Physical Units & Dimensional Scaling (100 Mpc, 32^3 grid, H0 = 70 km/s/Mpc)
        self.units = CosmologicalUnits(
            box_size_mpc=self.box_size_mpc,
            grid_resolution=self.grid_size,
            c_code=2.5,
            h0_km_s_mpc=70.0
        )

        # Physical Constants derived from CosmologicalUnits & Non-Equilibrium Thermodynamics
        self.DT = self.units.DT
        self.C_LIGHT = self.units.c_code
        self.KAPPA = self.units.get_cosmological_effective_kappa()
        self.H_0 = self.units.H_0
        self.G_CONST = self.units.G_CONST
        self.CS2 = self.units.CS2_BASE
        self.DIFFUSION_COEFF = self.units.DIFFUSION_BASE
        self.LANDAUER_DECAY = self.units.LANDAUER_BASE
        self.INFLATION_BOOST = self.units.INFLATION_BOOST
        self.ZETA_BEKENSTEIN = self.units.ZETA_BEKENSTEIN
        self.MASS_THRESHOLD = self.units.MASS_THRESHOLD
        self.M0_CORE = self.units.M0_CORE
        self.A_LOCAL_UNIVERSE = self.units.A_LOCAL_UNIVERSE
        self.A_MAX_CONFORMAL = self.units.A_MAX_CONFORMAL

        # Spatial Fourier Mesh on active backend (CPU/GPU)
        kx = 2.0 * np.pi * np.fft.fftfreq(self.grid_size)[:, None, None].astype(np.float32)
        ky = 2.0 * np.pi * np.fft.fftfreq(self.grid_size)[None, :, None].astype(np.float32)
        kz = 2.0 * np.pi * np.fft.fftfreq(self.grid_size)[None, None, :].astype(np.float32)
        k2_cpu = kx**2 + ky**2 + kz**2
        k2_cpu[0, 0, 0] = 1.0  # Regularize DC mode
        self.k2 = self.xp.asarray(k2_cpu)
        del kx, ky, kz, k2_cpu

        # High-Definition Mollweide Celestial Sphere Grid (90x180 = 16,200 samples)
        self.n_lat, self.n_lon = 90, 180
        self.lats = np.linspace(-np.pi / 2.0, np.pi / 2.0, self.n_lat)
        self.lons = np.linspace(-np.pi, np.pi, self.n_lon)
        self.LON, self.LAT = np.meshgrid(self.lons, self.lats)

        self.n_los_x = np.cos(self.LAT) * np.cos(self.LON)
        self.n_los_y = np.cos(self.LAT) * np.sin(self.LON)
        self.n_los_z = np.sin(self.LAT)

        self.r_obs = self.grid_size / 2.2
        self.cx_obs = self.grid_size / 2.0
        self.cy_obs = self.grid_size / 2.0
        self.cz_obs = self.grid_size / 2.0

        self.coords_cmb_x = (self.cx_obs + self.r_obs * self.n_los_x) % self.grid_size
        self.coords_cmb_y = (self.cy_obs + self.r_obs * self.n_los_y) % self.grid_size
        self.coords_cmb_z = (self.cz_obs + self.r_obs * self.n_los_z) % self.grid_size

        # Simulation State
        self.is_running = True
        self.steps_per_frame = self.initial_speed
        self.total_steps = 0
        self.t_coord = 0.0
        self.eon = 1
        self.scale_factor = 1.0
        self.s_bh_val = 0.0
        self.s_crit = 3500.0
        self.progress = 0.0
        self.mass_frac_val = 0.0
        self.last_bounce_step = 0
        self.eon_start_walltime = time.time()
        self.saved_epochs = set()

        # Initialize or Resume with grid-resolution tagged checkpoints
        latest_checkpoint = os.path.join(self.checkpoint_dir, f"latest_g{self.grid_size}.npz")
        if not os.path.exists(latest_checkpoint) and self.grid_size == 32:
            legacy_checkpoint = os.path.join(self.checkpoint_dir, "latest.npz")
            if os.path.exists(legacy_checkpoint):
                latest_checkpoint = legacy_checkpoint

        if force_reset:
            self.reset_simulation(archive_existing=True)
        elif auto_resume and os.path.exists(latest_checkpoint):
            self.load_checkpoint(latest_checkpoint)
        else:
            self._init_primordial_state()

    def to_cpu(self, arr):
        """Converts an array from GPU/CuPy or CPU memory to standard host NumPy ndarray."""
        if arr is None:
            return None
        if hasattr(arr, 'get'):
            return arr.get()
        return np.asarray(arr)

    @property
    def tau_physical(self) -> np.ndarray:
        """
        Physical emergent proper time along timelike trajectories:
        tau_phys(x, t) = t_coord + Delta_tau(x, t)
        Satisfies: d(tau_phys)/dt = 1 + kappa_0 * sigma_total(x, t).
        """
        return self.to_cpu(self.t_coord + self.tau).astype(np.float32)

    @property
    def tau_excess(self) -> np.ndarray:
        """
        Dissipative excess odometer tensor:
        Delta_tau(x, t) = int_0^t kappa * sigma_total(x, t') dt'.
        Vanishes in cosmic voids (sigma -> 0) and saturates inside virialized cores.
        """
        return self.to_cpu(self.tau).astype(np.float32)

    @property
    def d_tau_physical_dt(self) -> np.ndarray:
        """Instantaneous emergence rate of physical proper time: d(tau_phys)/dt = 1 + kappa * sigma."""
        return self.to_cpu(1.0 + self.d_tau_dt).astype(np.float32)

    @property
    def X(self) -> np.ndarray:
        """On-demand 3D X coordinate grid (CPU/GPU backend)."""
        return self.to_cpu(self.xp.arange(self.grid_size)[:, None, None] * self.xp.ones((1, self.grid_size, self.grid_size)))

    @property
    def Y(self) -> np.ndarray:
        """On-demand 3D Y coordinate grid (CPU/GPU backend)."""
        return self.to_cpu(self.xp.arange(self.grid_size)[None, :, None] * self.xp.ones((self.grid_size, 1, self.grid_size)))

    @property
    def Z(self) -> np.ndarray:
        """On-demand 3D Z coordinate grid (CPU/GPU backend)."""
        return self.to_cpu(self.xp.arange(self.grid_size)[None, None, :] * self.xp.ones((self.grid_size, self.grid_size, 1)))

    @property
    def p_k(self) -> np.ndarray:
        """On-demand CDM power spectrum array."""
        return self._compute_p_k()

    def _compute_p_k(self):
        """Computes Harrison-Zel'dovich + CDM turnover power spectrum P(k) on demand."""
        k_mod = self.xp.sqrt(self.k2)
        k_eq = 2.0 * np.pi * 3.5 / float(self.grid_size)
        q = k_mod / k_eq
        p_k = self.xp.zeros_like(self.k2)
        mask_k = k_mod > 0.0
        p_k[mask_k] = q[mask_k] / ((1.0 + (q[mask_k]**1.5))**2.0)
        return p_k

    def _init_primordial_state(self):
        """Initializes primordial cosmological fields with seeded Gaussian Random Field."""
        self.t_coord = 0.0
        # Primordial Gaussian Random Field (GRF) with scale-invariant Harrison-Zel'dovich power spectrum P(k)
        noise_raw = self.rng.standard_normal((self.grid_size, self.grid_size, self.grid_size)).astype(np.float32)
        noise_fft = self.xp.fft.fftn(self.xp.asarray(noise_raw))
        del noise_raw
        p_k = self._compute_p_k()
        fluct = self.xp.real(self.xp.fft.ifftn(noise_fft * self.xp.sqrt(p_k)))
        del p_k, noise_fft
        fluct_mean = float(self.to_cpu(self.xp.mean(fluct)))
        fluct_std = float(self.to_cpu(self.xp.std(fluct)))
        fluct = (fluct - fluct_mean) / max(1e-4, fluct_std) * 0.45

        # Organic primordial matter density field: rho_mean = 1.0 + delta_rho(k)
        self.rho = self.xp.maximum(0.05, 1.0 + fluct).astype(self.xp.float32)
        del fluct
        self.T = (12.0 * (self.rho**0.5) + 2.73).astype(self.xp.float32)
        self.I = self.xp.clip((self.rho - 0.5) / 2.5, 0.0, 1.0).astype(self.xp.float32)
        # self.tau stores the dissipative excess Delta tau (odometer)
        # Holographic Gravity Bootstrap (Eon 1 Synthetic Prior-Eon Fossil)
        # We assume Eon 0's thermal peaks highly correlated with Eon 1's primordial density.
        delta_rho_init = self.rho - 1.0
        HOLOGRAPHIC_BOOTSTRAP_FACTOR = 50.0
        synthetic_tau = (self.xp.maximum(0.0, delta_rho_init) * HOLOGRAPHIC_BOOTSTRAP_FACTOR).astype(self.xp.float32)
        
        self.tau = synthetic_tau.copy()
        self.tau_eon_start = synthetic_tau.copy()
        self.d_tau_dt = self.xp.zeros((self.grid_size, self.grid_size, self.grid_size), dtype=self.xp.float32)

        self.v_x = self.xp.zeros((self.grid_size, self.grid_size, self.grid_size), dtype=self.xp.float32)
        self.v_y = self.xp.zeros((self.grid_size, self.grid_size, self.grid_size), dtype=self.xp.float32)
        self.v_z = self.xp.zeros((self.grid_size, self.grid_size, self.grid_size), dtype=self.xp.float32)

    def reset_simulation(self, archive_existing=True):
        """
        Safely resets simulation back to primordial Eon 1.
        Archives previous run data to checkpoints/archive_<timestamp>/ for data safety.
        """
        if archive_existing and os.path.exists(self.checkpoint_dir):
            timestamp_str = time.strftime("%Y%m%d_%H%M%S")
            archive_path = os.path.join(self.checkpoint_dir, f"archive_{timestamp_str}")
            os.makedirs(archive_path, exist_ok=True)

            # Move older checkpoints and history to archive
            for item in os.listdir(self.checkpoint_dir):
                if item.startswith("archive_"):
                    continue
                src = os.path.join(self.checkpoint_dir, item)
                dst = os.path.join(archive_path, item)
                try:
                    if os.path.isdir(src):
                        shutil.copytree(src, dst, dirs_exist_ok=True)
                        shutil.rmtree(src)
                    else:
                        shutil.copy2(src, dst)
                        os.remove(src)
                except Exception:
                    pass

        # Re-create clean snapshots dir & empty history
        os.makedirs(self.snapshots_dir, exist_ok=True)
        with open(self.history_file, 'w', encoding='utf-8') as f:
            json.dump([], f)

        # Reset all cosmological variables to primordial Eon 1
        self.eon = 1
        self.scale_factor = 1.0
        self.total_steps = 0
        self.s_bh_val = 0.0
        self.progress = 0.0
        self.mass_frac_val = 0.0
        self.last_bounce_step = 0
        self.eon_start_walltime = time.time()

        self._init_primordial_state()
        self.save_checkpoint()
        return True

    def _sample_sphere_trilinear(self, arr, cx, cy, cz):
        """Continuous trilinear interpolation across periodic 3-torus boundary conditions."""
        cx_xp = self.xp.asarray(cx)
        cy_xp = self.xp.asarray(cy)
        cz_xp = self.xp.asarray(cz)

        x0 = self.xp.floor(cx_xp).astype(int) % self.grid_size
        x1 = (x0 + 1) % self.grid_size
        y0 = self.xp.floor(cy_xp).astype(int) % self.grid_size
        y1 = (y0 + 1) % self.grid_size
        z0 = self.xp.floor(cz_xp).astype(int) % self.grid_size
        z1 = (z0 + 1) % self.grid_size

        xd = cx_xp - self.xp.floor(cx_xp)
        yd = cy_xp - self.xp.floor(cy_xp)
        zd = cz_xp - self.xp.floor(cz_xp)

        c00 = arr[x0, y0, z0] * (1.0 - xd) + arr[x1, y0, z0] * xd
        c01 = arr[x0, y0, z1] * (1.0 - xd) + arr[x1, y0, z1] * xd
        c10 = arr[x0, y1, z0] * (1.0 - xd) + arr[x1, y1, z0] * xd
        c11 = arr[x0, y1, z1] * (1.0 - xd) + arr[x1, y1, z1] * xd

        c0 = c00 * (1.0 - yd) + c10 * yd
        c1 = c01 * (1.0 - yd) + c11 * yd

        return c0 * (1.0 - zd) + c1 * zd

    def _generate_phase_locked_fluctuations(self, tau: np.ndarray, alpha_mem: float = 0.35) -> np.ndarray:
        """
        Generates primordial matter fluctuations via Holographic Phase-Locking in Fourier Space.
        Combines the phase angle of the multi-eon fossil tensor tau(x) with quantum vacuum fluctuations.
        """
        tau_xp = self.xp.asarray(tau, dtype=self.xp.float32)
        tau_fft = self.xp.fft.fftn(tau_xp)
        theta_fossil = self.xp.angle(tau_fft)
        del tau_fft

        noise_raw = self.rng.standard_normal((self.grid_size, self.grid_size, self.grid_size)).astype(np.float32)
        noise_fft = self.xp.fft.fftn(self.xp.asarray(noise_raw))
        del noise_raw
        theta_quantum = self.xp.angle(noise_fft)
        del noise_fft

        z_fossil = self.xp.exp(1j * theta_fossil)
        del theta_fossil
        z_quantum = self.xp.exp(1j * theta_quantum)
        del theta_quantum
        z_bounce = alpha_mem * z_fossil + (1.0 - alpha_mem) * z_quantum
        del z_fossil, z_quantum
        theta_bounce = self.xp.angle(z_bounce)
        del z_bounce

        p_k = self._compute_p_k()
        synthesized_fft = self.xp.sqrt(p_k) * self.xp.exp(1j * theta_bounce)
        del p_k, theta_bounce
        synthesized_fft[0, 0, 0] = 0.0 + 0.0j

        fluct = self.xp.real(self.xp.fft.ifftn(synthesized_fft))
        del synthesized_fft
        std_val = float(self.to_cpu(self.xp.std(fluct)))
        mean_val = float(self.to_cpu(self.xp.mean(fluct)))
        if std_val > 1e-6:
            fluct = (fluct - mean_val) / std_val * 0.35

        return fluct.astype(self.xp.float32)

    def _trigger_white_hole_eon_3d(self):
        """
        Detonates a localized white hole quantum bounce (Carlo Rovelli / Planck Star metric inversion)
        when a singular gravitational core exceeds the critical mass fraction threshold (Route A).
        """
        if float(self.to_cpu(self.xp.max(self.rho))) > 1.2:
            flat_idx = int(self.to_cpu(self.xp.argmax(self.rho)))
        else:
            flat_idx = int(self.to_cpu(self.xp.argmax(self.tau)))

        x0, y0, z0 = np.unravel_index(flat_idx, (self.grid_size, self.grid_size, self.grid_size))

        dx = ((self.xp.arange(self.grid_size) - x0 + self.grid_size / 2.0) % self.grid_size - self.grid_size / 2.0)[:, None, None]
        dy = ((self.xp.arange(self.grid_size) - y0 + self.grid_size / 2.0) % self.grid_size - self.grid_size / 2.0)[None, :, None]
        dz = ((self.xp.arange(self.grid_size) - z0 + self.grid_size / 2.0) % self.grid_size - self.grid_size / 2.0)[None, None, :]
        r = self.xp.sqrt(dx**2 + dy**2 + dz**2)
        r_safe = self.xp.maximum(0.8, r)

        fluct_new = self._generate_phase_locked_fluctuations(self.tau, alpha_mem=0.35)

        primordial_blast = 3.5 * self.xp.exp(-r**2 / 16.0)
        thermal_reheating = 85.0 * self.xp.exp(-r**2 / 20.0) + 2.73

        v_exp_mag = 2.4 * self.xp.exp(-r**2 / 25.0) * (r / self.grid_size)
        v_x_new = v_exp_mag * (dx / r_safe)
        v_y_new = v_exp_mag * (dy / r_safe)
        v_z_new = v_exp_mag * (dz / r_safe)
        del dx, dy, dz, r, r_safe, v_exp_mag

        rho_new = self.xp.clip(1.0 + fluct_new + primordial_blast, 0.05, 12.0)
        T_new = self.xp.clip(thermal_reheating + 15.0 * self.xp.abs(fluct_new), 2.73, 2000.0)
        del fluct_new, primordial_blast, thermal_reheating

        return rho_new, v_x_new, v_y_new, v_z_new, T_new

    def _trigger_conformal_eon_3d(self):
        """
        Executes Roger Penrose's Conformal Cyclic Cosmology (CCC) transition across the spacelike hypersurface I+.
        Resets the universe globally and homogeneously across the entire 3D lattice, imprinting the fossil
        holographic memory of the prior eon's proper time tensor tau into primordial multi-scale perturbations.
        """
        # 1. Synthesize global scale-invariant density fluctuations with holographic phase-locking
        fluct_new = self._generate_phase_locked_fluctuations(self.tau, alpha_mem=0.35)

        # 2. Global homogeneous matter density field (centered at mean rho = 1.0)
        rho_new = self.xp.clip(1.0 + fluct_new, 0.05, 5.0)

        # 3. Global isotropic CMB reheating temperature (base 2.73 K + Sachs-Wolfe plasma fluctuations)
        T_new = self.xp.clip(2.73 + 12.0 * self.xp.abs(fluct_new), 2.73, 500.0)

        # 4. Primordial velocity fields via Zel'dovich gravitational perturbation equation: v = -grad(Phi_primordial) * dt
        delta_rho = rho_new - self.xp.mean(rho_new)
        phi_prim_fft = self.xp.fft.fftn(delta_rho)
        del delta_rho
        phi_prim_fft *= (-4.0 * np.pi * self.G_CONST) / (self.k2 * 1.0)
        phi_prim_fft[0, 0, 0] = 0.0
        phi_prim = self.xp.real(self.xp.fft.ifftn(phi_prim_fft))
        del phi_prim_fft

        # Zel'dovich velocity: v = -grad(Phi) * t_primordial
        v_x_new = -0.35 * self._grad_axis(phi_prim, 0)
        v_y_new = -0.35 * self._grad_axis(phi_prim, 1)
        v_z_new = -0.35 * self._grad_axis(phi_prim, 2)
        del phi_prim, fluct_new

        return rho_new.astype(self.xp.float32), v_x_new.astype(self.xp.float32), v_y_new.astype(self.xp.float32), v_z_new.astype(self.xp.float32), T_new.astype(self.xp.float32)

    def _grad_axis(self, arr, axis):
        """Computes central spatial difference along a single axis."""
        return 0.5 * (self.xp.roll(arr, -1, axis=axis) - self.xp.roll(arr, 1, axis=axis))

    def _compute_rhs(self, a_factor, rho, vx, vy, vz, T, I, tau):
        xp = self.xp
        # 1. H_eff
        if a_factor < 1.05:
            H_eff = self.H_0 * (1.0 + self.INFLATION_BOOST * np.exp(-(a_factor - 1.0) / 0.015))
        else:
            H_eff = self.H_0

        # 2. Phi (Holographic Gravity Boost)
        # tau_eon_start acts as localized Apparent Dark Matter, overcoming Jeans mass without real dark particles
        HOLOGRAPHIC_COUPLING = 1.0
        delta_rho = rho - xp.mean(rho)
        delta_rho_eff = delta_rho + HOLOGRAPHIC_COUPLING * self.tau_eon_start
        
        phi_fft = xp.fft.fftn(delta_rho_eff)
        phi_fft *= (-4.0 * np.pi * self.G_CONST) / (self.k2 * max(1.0, float(a_factor)))
        phi_fft[0, 0, 0] = 0.0
        phi = xp.real(xp.fft.ifftn(phi_fft))
        
        # 3. Hydro
        cs2 = self.units.compute_sound_speed_sq(T, base_cs2=self.CS2, xp=xp)
        hubble_damp = 2.0 * H_eff / max(1.0, float(a_factor))
        
        div_rho = xp.zeros_like(rho)
        div_px = xp.zeros_like(rho)
        div_py = xp.zeros_like(rho)
        div_pz = xp.zeros_like(rho)
        
        c_num = 0.04
        grad_phi_sq = xp.zeros_like(rho)
        
        px, py, pz = rho * vx, rho * vy, rho * vz
        
        for axis, v_ax in enumerate([vx, vy, vz]):
            # Mass flux
            f_rho = rho * v_ax
            r_f_p = xp.roll(f_rho, -1, axis=axis)
            r_r_p = xp.roll(rho, -1, axis=axis)
            fp_rho = 0.5 * (f_rho + r_f_p) - 0.5 * c_num * (r_r_p - rho)
            
            r_f_m = xp.roll(f_rho, 1, axis=axis)
            r_r_m = xp.roll(rho, 1, axis=axis)
            fm_rho = 0.5 * (r_f_m + f_rho) - 0.5 * c_num * (rho - r_r_m)
            div_rho += (fp_rho - fm_rho)
            
            # Mom fluxes
            for target_p, p_ax in zip([div_px, div_py, div_pz], [px, py, pz]):
                f_p_ax = p_ax * v_ax
                r_fp_p = xp.roll(f_p_ax, -1, axis=axis)
                r_rp_p = xp.roll(p_ax, -1, axis=axis)
                fp_pax = 0.5 * (f_p_ax + r_fp_p) - 0.5 * c_num * (r_rp_p - p_ax)
                
                r_fm_p = xp.roll(f_p_ax, 1, axis=axis)
                r_rm_p = xp.roll(p_ax, 1, axis=axis)
                fm_pax = 0.5 * (r_fm_p + f_p_ax) - 0.5 * c_num * (p_ax - r_rm_p)
                target_p += (fp_pax - fm_pax)
                
            g_phi = self._grad_axis(phi, axis)
            grad_phi_sq += g_phi**2
            
            # Source term for momentum: -grad(P) - rho * grad(Phi)
            g_P = self._grad_axis(cs2 * rho, axis)
            source_p = -g_P - rho * g_phi - hubble_damp * (px if axis==0 else (py if axis==1 else pz))
            
            if axis == 0: div_px = div_px - source_p
            elif axis == 1: div_py = div_py - source_p
            else: div_pz = div_pz - source_p

        d_rho_dt = -div_rho
        d_px_dt = -div_px
        d_py_dt = -div_py
        d_pz_dt = -div_pz
        
        # 4. Thermal
        kappa = self.units.compute_spitzer_conductivity(T, rho, base_k=self.DIFFUSION_COEFF, xp=xp)
        laplacian_T = (xp.roll(T, 1, 0) + xp.roll(T, -1, 0) +
                       xp.roll(T, 1, 1) + xp.roll(T, -1, 1) +
                       xp.roll(T, 1, 2) + xp.roll(T, -1, 2) - 6.0 * T)
                       
        T_ad = 2.73 + 12.0 * (rho ** 0.67) / max(1.0, float(a_factor)**0.5)
        comp_heating = 0.05 * xp.maximum(0.0, d_rho_dt) * T
        hubble_cooling = (H_eff / max(1.0, float(a_factor))) * (T - 2.73)
        d_T_dt = 0.25 * kappa * laplacian_T + comp_heating - hubble_cooling + 0.1 * (T_ad - T)
        
        # 5. Onsager Entropy & Tau (with Sub-Grid Virial Shock Enhancer)
        # Coarse grids mathematically smear sub-megaparsec thermal shocks.
        # We amplify local entropy production in non-linear collapsed halos (overdensities)
        # to correctly predict asymptotic macroscopic time-dilation.
        inv_T = 1.0 / T
        sigma_th = xp.zeros_like(T)
        for axis in range(3):
            sigma_th += (-kappa * self._grad_axis(T, axis)) * self._grad_axis(inv_T, axis)
        sigma_th = xp.maximum(0.0, sigma_th)
        sigma_grav = (rho * grad_phi_sq) / (T * 50.0)
        
        mean_rho = xp.maximum(xp.mean(rho), 1e-10)
        overdensity = rho / mean_rho
        
        # Exponential sub-grid scaling triggered inside virialized clusters.
        # The critical overdensity threshold scales inversely with voxel volume.
        # Calibrated for DX_ref ~ 3.9 Mpc (Grid 128 in 500 Mpc box).
        ref_dx = 500.0 / 128.0
        current_dx = float(self.box_size_mpc) / float(self.grid_size)
        delta_threshold = max(2.0, 50.0 * (ref_dx / current_dx)**3)
        
        # Subgrid Virial Enhancer (Clipped to prevent numerical runaway singularities)
        subgrid_factor = 1.0 + xp.clip(xp.maximum(0.0, (overdensity / delta_threshold)**2), 0.0, 1000.0) * 100.0
        
        sigma_total = (sigma_th + sigma_grav) * subgrid_factor
        d_tau_dt = self.KAPPA * sigma_total
        
        # 6. Information
        div_I = xp.zeros_like(I)
        for axis, v_ax in enumerate([vx, vy, vz]):
            div_I += self._grad_axis(I * v_ax, axis)
        lap_I = (xp.roll(I, 1, 0) + xp.roll(I, -1, 0) +
                 xp.roll(I, 1, 1) + xp.roll(I, -1, 1) +
                 xp.roll(I, 1, 2) + xp.roll(I, -1, 2) - 6.0 * I)
        sustenance = 0.6 * sigma_total * (rho / xp.mean(rho))
        decay = self.units.compute_landauer_decay(T, base_decay=self.LANDAUER_DECAY, xp=xp)
        d_I_dt = -div_I + 0.02 * lap_I + (sustenance - decay * I)
        
        return d_rho_dt, d_px_dt, d_py_dt, d_pz_dt, d_T_dt, d_I_dt, d_tau_dt, H_eff

    def step(self):
        """Executes a single SSP-RK2 cosmological differential step."""
        xp = self.xp
        dt = self.DT
        
        # Current state U^n
        a_n = self.scale_factor
        rho_n = self.rho
        vx_n, vy_n, vz_n = self.v_x, self.v_y, self.v_z
        T_n = self.T
        I_n = self.I
        tau_n = self.tau
        
        # Stage 1: U^{(1)} = U^n + dt * L(U^n)
        drho, dpx, dpy, dpz, dT, dI, dtau, H_eff = self._compute_rhs(a_n, rho_n, vx_n, vy_n, vz_n, T_n, I_n, tau_n)
        
        a_1 = a_n + dt * H_eff
        
        # Mass clipping diagnostic
        rho_candidate = rho_n + dt * drho
        mass_before = float(self.to_cpu(xp.sum(rho_candidate)))
        rho_1 = xp.clip(rho_candidate, 1e-4, 1e5)
        mass_after = float(self.to_cpu(xp.sum(rho_1)))
        self.mass_clip_error = abs(mass_after - mass_before)
        
        # Limit momentum via causal limit (v <= c)
        px_1 = rho_n * vx_n + dt * dpx
        py_1 = rho_n * vy_n + dt * dpy
        pz_1 = rho_n * vz_n + dt * dpz
        
        vx_1 = px_1 / rho_1
        vy_1 = py_1 / rho_1
        vz_1 = pz_1 / rho_1
        
        v_mag_1 = xp.sqrt(vx_1**2 + vy_1**2 + vz_1**2)
        v_limit_1 = xp.maximum(1.0, v_mag_1 / self.C_LIGHT)
        vx_1 /= v_limit_1; vy_1 /= v_limit_1; vz_1 /= v_limit_1
        
        T_1 = xp.clip(T_n + dt * dT, 2.73, 1e8)
        I_1 = xp.clip(I_n + dt * dI, 0.0, 1.0)
        tau_1 = tau_n + dt * dtau
        
        # Stage 2: U^{n+1} = 0.5 * U^n + 0.5 * U^{(1)} + 0.5 * dt * L(U^{(1)})
        drho2, dpx2, dpy2, dpz2, dT2, dI2, dtau2, H_eff2 = self._compute_rhs(a_1, rho_1, vx_1, vy_1, vz_1, T_1, I_1, tau_1)
        
        self.scale_factor = 0.5 * a_n + 0.5 * a_1 + 0.5 * dt * H_eff2
        self.rho = xp.clip(0.5 * rho_n + 0.5 * rho_1 + 0.5 * dt * drho2, 1e-4, 1e5)
        
        px_2 = 0.5 * (rho_n * vx_n) + 0.5 * (rho_1 * vx_1) + 0.5 * dt * dpx2
        py_2 = 0.5 * (rho_n * vy_n) + 0.5 * (rho_1 * vy_1) + 0.5 * dt * dpy2
        pz_2 = 0.5 * (rho_n * vz_n) + 0.5 * (rho_1 * vz_1) + 0.5 * dt * dpz2
        
        self.v_x = px_2 / self.rho
        self.v_y = py_2 / self.rho
        self.v_z = pz_2 / self.rho
        
        v_mag_2 = xp.sqrt(self.v_x**2 + self.v_y**2 + self.v_z**2)
        v_limit_2 = xp.maximum(1.0, v_mag_2 / self.C_LIGHT)
        self.v_x /= v_limit_2; self.v_y /= v_limit_2; self.v_z /= v_limit_2
        
        self.T = xp.clip(0.5 * T_n + 0.5 * T_1 + 0.5 * dt * dT2, 2.73, 1e8)
        self.I = xp.clip(0.5 * I_n + 0.5 * I_1 + 0.5 * dt * dI2, 0.0, 1.0)
        
        # Bekenstein-Hawking Saturation Limit for Virialized Cores
        self.tau = xp.clip(0.5 * tau_n + 0.5 * tau_1 + 0.5 * dt * dtau2, 0.0, 1e7)
        self.d_tau_dt = xp.clip(dtau2, -1e6, 1e7)
        
        self.t_coord += dt

        # 8. Bekenstein Quantum Saturation & Virialized Structure Diagnostics
        tau_current_eon = self.tau - self.tau_eon_start
        total_mass = float(self.to_cpu(self.xp.sum(self.rho)))
        # Virialized compact core criterion (overdensity delta_rho / rho_bar >= 2.0 -> rho >= 3.0)
        core_mask = self.rho > 3.0
        core_mass = float(self.to_cpu(self.xp.sum(self.rho[core_mask])))
        self.mass_frac_val = core_mass / max(1.0, total_mass)

        self.s_bh_val = float(self.to_cpu(self.xp.max(tau_current_eon)))
        self.s_crit = self.units.compute_bekenstein_entropy_limit(
            mass_core=core_mass,
            m0_ref=self.M0_CORE,
            zeta_base=self.ZETA_BEKENSTEIN
        )

        # Global Conformal Boundary Progress: Normalized expansion towards Penrose asymptotic heat death (a -> 7.0)
        self.progress = max(0.0, min(1.0, (self.scale_factor - 1.0) / (self.A_MAX_CONFORMAL - 1.0)))
        self.total_steps += 1

        # Check and save scientific epoch checkpoints across cosmological history
        g_tag = f"_g{self.grid_size}"
        if self.scale_factor >= 1.00 and "cmb" not in self.saved_epochs:
            self.saved_epochs.add("cmb")
            self.save_checkpoint(os.path.join(self.checkpoint_dir, f"cmb_eon_{self.eon}{g_tag}.npz"))
        if self.scale_factor >= 1.50 and "dawn" not in self.saved_epochs:
            self.saved_epochs.add("dawn")
            self.save_checkpoint(os.path.join(self.checkpoint_dir, f"dawn_eon_{self.eon}{g_tag}.npz"))
        if self.scale_factor >= 2.00 and "bao" not in self.saved_epochs:
            self.saved_epochs.add("bao")
            self.save_checkpoint(os.path.join(self.checkpoint_dir, f"bao_eon_{self.eon}{g_tag}.npz"))
        if self.scale_factor >= 3.00 and "clusters" not in self.saved_epochs:
            self.saved_epochs.add("clusters")
            self.save_checkpoint(os.path.join(self.checkpoint_dir, f"clusters_eon_{self.eon}{g_tag}.npz"))
        if self.scale_factor >= 4.50 and "pantheon" not in self.saved_epochs:
            self.saved_epochs.add("pantheon")
            self.save_checkpoint(os.path.join(self.checkpoint_dir, f"pantheon_eon_{self.eon}{g_tag}.npz"))

        # Cosmological Phase Diagram: Global CCC Saturation
        is_conformal_bounce = (self.scale_factor >= self.A_MAX_CONFORMAL)

        if is_conformal_bounce:
            self._handle_bounce(transition_type="CCC (Muerte Térmica)")

    def _handle_bounce(self, transition_type="CCC (Muerte Térmica)"):
        """Processes transition to next eon, logs history, archives full snapshot, and notifies Telegram."""
        eon_duration_wall = time.time() - self.eon_start_walltime
        eon_steps = self.total_steps - self.last_bounce_step

        # Save Full Visual Snapshot of the Completed Eon at Bounce
        final_snapshot = self.get_visual_payload()
        final_snapshot["snapshot_meta"] = {
            "id": f"eon_{self.eon}_g{self.grid_size}",
            "label": f"📷 Fin Eón {self.eon} ({self.grid_size}³) [{transition_type}]",
            "type": "eon_bounce",
            "transition": transition_type,
            "eon": self.eon,
            "grid_size": self.grid_size,
            "scale_factor": round(float(self.scale_factor), 3),
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        snapshot_path = os.path.join(self.snapshots_dir, f"snapshot_eon_{self.eon}_g{self.grid_size}.json")
        try:
            with open(snapshot_path, "w", encoding="utf-8") as f:
                json.dump(final_snapshot, f)
        except Exception:
            pass

        # Locate central attractor
        flat_idx_rho = int(self.to_cpu(self.xp.argmax(self.rho)))
        bx, by, bz = np.unravel_index(flat_idx_rho, self.rho.shape)

        # Generate Observational Reports (CMB Mollweide Map + Planck Power Spectrum)
        observational_metrics = {}
        h0_pred = self.units.H0_PLANCK_BASELINE
        c2_val = 0.0
        c3_val = 0.0
        c2_c3_ratio = 1.0
        planck_chi2_val = 0.0
        try:
            from experiments.compare_planck import generate_eon_observational_report
            observational_metrics = generate_eon_observational_report(self, output_dir=self.snapshots_dir)
            h0_pred = float(observational_metrics.get("h0_predicted", self.units.H0_PLANCK_BASELINE))
            c2_val = float(observational_metrics.get("quadrupole_C2", 0.0))
            c3_val = float(observational_metrics.get("octopole_C3", 0.0))
            c2_c3_ratio = float(observational_metrics.get("ratio_C2_C3", 1.0))
            planck_chi2_val = float(observational_metrics.get("planck_chi2", 0.0))
        except Exception as ex:
            observational_metrics = {"error": str(ex)}

        # Record Enriched History Entry
        history_entry = {
            "eon": self.eon,
            "transition": transition_type,
            "final_scale_factor": round(float(self.scale_factor), 3),
            "peak_s_bh": round(float(self.s_bh_val), 1),
            "s_crit": round(float(self.s_crit), 1),
            "core_mass_fraction": round(float(self.mass_frac_val * 100.0), 2),
            "fossil_odometer_total": round(float(self.to_cpu(self.xp.max(self.tau))), 1),
            "h0_predicted": round(h0_pred, 2),
            "quadrupole_C2": round(c2_val, 4),
            "octopole_C3": round(c3_val, 4),
            "ratio_C2_C3": round(c2_c3_ratio, 3),
            "planck_chi2": round(planck_chi2_val, 2),
            "attractor": {"x": int(bx), "y": int(by), "z": int(bz)},
            "eon_steps": eon_steps,
            "walltime_seconds": round(eon_duration_wall, 1),
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "observational": observational_metrics
        }
        self._append_history(history_entry)

        # Trigger Telegram Alert if configured
        self.notifier.check_and_notify_eon(history_entry)

        # Save Completed Eon Final Archive (Epoch 6)
        g_tag = f"_g{self.grid_size}"
        self.save_checkpoint(os.path.join(self.checkpoint_dir, f"eon_{self.eon}{g_tag}.npz"))

        # Increment Eon
        self.eon += 1
        self.scale_factor = 1.0
        self.saved_epochs = set()
        self.tau_eon_start = self.tau.copy()
        if transition_type == "Rebote Gravitatorio (Túnel Cuántico)":
            self.rho, self.v_x, self.v_y, self.v_z, self.T = self._trigger_white_hole_eon_3d()
        else:
            self.rho, self.v_x, self.v_y, self.v_z, self.T = self._trigger_conformal_eon_3d()
        self.I = self.xp.clip((self.rho - 0.5) / 2.5, 0.0, 1.0)
        self.progress = 0.0
        self.last_bounce_step = self.total_steps
        self.eon_start_walltime = time.time()

        # Save Checkpoints
        self.save_checkpoint()

    def step_batch(self, n_steps=None):
        """Executes a batch of steps."""
        steps = n_steps if n_steps is not None else self.steps_per_frame
        for _ in range(steps):
            self.step()

    def get_telemetry(self):
        """Returns structured real-time telemetry dictionary."""
        rho_cpu = self.to_cpu(self.rho)
        tau_cpu = self.to_cpu(self.tau)
        T_cpu = self.to_cpu(self.T)

        bx, by, bz = np.unravel_index(np.argmax(rho_cpu), rho_cpu.shape)
        z_slice = int(np.clip(bz, 0, self.grid_size - 1))

        if self.scale_factor < 1.05:
            era_str = "Fase de Inflación Cuántica Primordial"
        elif self.scale_factor < 2.5:
            era_str = "Era de Filamentos y Panqueques 3D"
        elif self.scale_factor < 4.5:
            era_str = "Era de Cúmulos Virializados y Galaxias 3D"
        elif self.scale_factor < 6.0:
            era_str = "Era de Expansión Acelerada y Desacople"
        else:
            era_str = "Fase Asintótica Pre-Rebote 3D (Límite Conforme CCC)"

        redshift = max(0.0, (self.A_LOCAL_UNIVERSE / max(0.01, float(self.scale_factor))) - 1.0)
        temp_norm = float(np.percentile(T_cpu, 99))
        temp_astro = temp_norm * 120.0

        if self.progress >= 0.95:
            status_banner = "Transición Conforme CCC Inminente (Nuevo Eón)"
        elif self.scale_factor >= 5.5:
            status_banner = "Dilución Asintótica (Fase Conforme de Penrose)"
        else:
            status_banner = "Evolución Hidrodinámica 3D"

        prog_label = f"Frontera Conforme CCC Eón {self.eon}"
        active_route = "conformal"
        p_conformal = float(self.progress)
        p_grav = float(min(1.0, self.s_bh_val / max(1.0, self.s_crit)))

        tau_phys_max = float(self.t_coord + np.max(tau_cpu))
        time_myr = float(self.units.time_code_to_myr(tau_phys_max))

        return {
            "eon": self.eon,
            "era": era_str,
            "scale_factor": round(float(self.scale_factor), 3),
            "redshift": round(float(redshift), 2),
            "temp_norm": round(temp_norm, 1),
            "temp_astro": round(temp_astro, 0),
            "c_light": round(self.C_LIGHT, 2),
            "mass_fraction": round(float(self.mass_frac_val * 100.0), 1),
            "s_bh": round(float(self.s_bh_val), 0),
            "s_crit": round(float(self.s_crit), 0),
            "tunnel_progress": round(float(self.progress * 100.0), 1),
            "progress_label": prog_label,
            "active_route": active_route,
            "p_grav": round(float(p_grav * 100.0), 1),
            "p_conformal": round(float(p_conformal * 100.0), 1),
            "fossil_odometer": round(tau_phys_max, 1),
            "time_myr": round(time_myr, 1),
            "grid_size": self.grid_size,
            "grid_voxels": self.grid_size ** 3,
            "box_size_mpc": self.units.box_size_mpc,
            "cell_size_mpc": round(self.units.cell_size_mpc, 3),
            "h0_kms_mpc": self.units.h0_kms_mpc,
            "kappa_0_planck": f"{self.units.get_fundamental_kappa_0():.3e}",
            "attractor": {"x": int(bx), "y": int(by), "z": int(bz)},
            "z_slice": z_slice,
            "total_steps": self.total_steps,
            "is_running": self.is_running,
            "steps_per_frame": self.steps_per_frame,
            "use_gpu": self.use_gpu,
            "hardware": "GPU (CuPy / CUDA)" if self.use_gpu else "CPU (NumPy / OpenBLAS)",
            "state_status": status_banner
        }

    def get_visual_payload(self):
        """Constructs data arrays for the 9-panel web dashboard."""
        rho_cpu = self.to_cpu(self.rho)
        tau_cpu = self.to_cpu(self.tau)
        T_cpu = self.to_cpu(self.T)
        d_tau_dt_cpu = self.to_cpu(self.d_tau_dt)
        I_cpu = self.to_cpu(self.I)

        bx, by, bz = np.unravel_index(np.argmax(rho_cpu), rho_cpu.shape)
        z_slice = int(np.clip(bz, 0, self.grid_size - 1))

        # 1. 3D Cosmic Web point cloud (dense filaments and cluster nodes)
        mean_r = float(np.mean(rho_cpu))
        p90 = float(np.percentile(rho_cpu, 90))
        threshold_rho = max(mean_r * 1.05, p90)
        xs, ys, zs = np.where(rho_cpu >= threshold_rho)
        points_3d = []
        if len(xs) > 0:
            step_stride = max(1, len(xs) // 350)
            for i in range(0, len(xs), step_stride):
                points_3d.append([
                    int(xs[i]), int(ys[i]), int(zs[i]),
                    round(float(rho_cpu[xs[i], ys[i], zs[i]]), 2)
                ])

        # 2. HD Mollweide CMB (90x180) with Physical Sachs-Wolfe, Plasma Perturbations & Doppler Shift
        tau_s = self.to_cpu(self._sample_sphere_trilinear(self.tau, self.coords_cmb_x, self.coords_cmb_y, self.coords_cmb_z))
        T_s = self.to_cpu(self._sample_sphere_trilinear(self.T, self.coords_cmb_x, self.coords_cmb_y, self.coords_cmb_z))
        rho_s = self.to_cpu(self._sample_sphere_trilinear(self.rho, self.coords_cmb_x, self.coords_cmb_y, self.coords_cmb_z))
        vx_s = self.to_cpu(self._sample_sphere_trilinear(self.v_x, self.coords_cmb_x, self.coords_cmb_y, self.coords_cmb_z))
        vy_s = self.to_cpu(self._sample_sphere_trilinear(self.v_y, self.coords_cmb_x, self.coords_cmb_y, self.coords_cmb_z))
        vz_s = self.to_cpu(self._sample_sphere_trilinear(self.v_z, self.coords_cmb_x, self.coords_cmb_y, self.coords_cmb_z))

        v_los = (vx_s * self.n_los_x + vy_s * self.n_los_y + vz_s * self.n_los_z) / self.C_LIGHT
        mean_T = max(1e-4, float(np.mean(T_cpu)))
        mean_rho = max(1e-4, float(np.mean(rho_cpu)))
        mean_tau = max(1e-4, float(np.mean(tau_cpu)) + 1.0)

        delta_T_int = (T_s - mean_T) / mean_T
        delta_rho = (rho_s - mean_rho) / mean_rho
        delta_tau = (tau_s - float(np.mean(tau_cpu))) / mean_tau

        cmb_raw = delta_T_int + (1.0 / 3.0) * delta_rho + v_los + 0.35 * delta_tau
        cmb_std = max(1e-4, float(np.std(cmb_raw)))
        cmb_norm = np.clip((cmb_raw - float(np.mean(cmb_raw))) / cmb_std, -2.5, 2.5)

        # 3. 2D Cross Sections (32x32) at z_slice
        tau_phys_slice = tau_cpu[:, :, z_slice] + self.t_coord
        slice_rho = np.round(rho_cpu[:, :, z_slice], 3).tolist()
        slice_rate = np.round(d_tau_dt_cpu[:, :, z_slice], 3).tolist()
        slice_index = np.round(I_cpu[:, :, z_slice], 3).tolist()
        slice_tau = np.round(tau_phys_slice, 1).tolist()
        slice_log_tau = np.round(np.log10(1.0 + np.maximum(0.0, tau_phys_slice)), 3).tolist()
        slice_temp = np.round(T_cpu[:, :, z_slice], 2).tolist()

        return {
            "telemetry": self.get_telemetry(),
            "points_3d": points_3d,
            "cmb": np.round(cmb_norm, 2).tolist(),
            "slice_rho": slice_rho,
            "slice_rate": slice_rate,
            "slice_index": slice_index,
            "slice_tau": slice_tau,
            "slice_log_tau": slice_log_tau,
            "slice_temp": slice_temp
        }

    def save_manual_snapshot(self):
        """Saves current state and visual payload as an interactive manual checkpoint snapshot."""
        # Save physical binary arrays
        self.save_checkpoint()

        # Save visual snapshot JSON with metadata
        timestamp_id = time.strftime("%Y%m%d_%H%M%S")
        time_display = time.strftime("%H:%M:%S")
        snapshot_id = f"manual_{timestamp_id}_g{self.grid_size}"
        
        payload = self.get_visual_payload()
        payload["snapshot_meta"] = {
            "id": snapshot_id,
            "label": f"💾 Guardado: Eón {self.eon} (a={self.scale_factor:.3f}, {time_display}) [{self.grid_size}³]",
            "type": "manual",
            "eon": self.eon,
            "grid_size": self.grid_size,
            "scale_factor": round(float(self.scale_factor), 3),
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }

        snapshot_path = os.path.join(self.snapshots_dir, f"snapshot_{snapshot_id}.json")
        try:
            with open(snapshot_path, "w", encoding="utf-8") as f:
                json.dump(payload, f)
        except Exception:
            pass

        return payload["snapshot_meta"]

    def get_snapshot(self, snapshot_id):
        """Retrieves an archived visual snapshot by ID or eon number."""
        str_id = str(snapshot_id).strip()
        
        # Candidate file names (matching both generic, tagged, and grid-suffixed)
        candidates = [
            f"snapshot_{str_id}.json",
            f"snapshot_{str_id}_g{self.grid_size}.json",
            f"snapshot_eon_{str_id}.json",
            f"snapshot_eon_{str_id}_g{self.grid_size}.json",
            f"snapshot_manual_{str_id}.json",
            f"snapshot_manual_{str_id}_g{self.grid_size}.json",
            f"{str_id}.json"
        ]

        for cand in candidates:
            path = os.path.join(self.snapshots_dir, cand)
            if os.path.exists(path):
                try:
                    with open(path, "r", encoding="utf-8") as f:
                        return json.load(f)
                except Exception:
                    pass
        return None

    @staticmethod
    def normalize_transition_name(transition_str):
        """Normalizes legacy and modern transition names into standard nomenclature."""
        if not transition_str:
            return "Rebote Gravitatorio (Túnel Cuántico)"
        t = str(transition_str).strip()
        if "CCC" in t or "Muerte" in t or "Conforme" in t:
            return "CCC (Muerte Térmica)"
        return "Rebote Gravitatorio (Túnel Cuántico)"

    def get_available_snapshots(self):
        """Returns list of all available snapshots with strictly 1 entry per eon and normalized transition labels."""
        if not os.path.exists(self.snapshots_dir):
            return []
        files = os.listdir(self.snapshots_dir)
        snapshots = []
        seen_eons = set()
        seen_manual_ids = set()

        for f in files:
            if f.startswith("snapshot_") and f.endswith(".json"):
                path = os.path.join(self.snapshots_dir, f)
                try:
                    with open(path, "r", encoding="utf-8") as sfile:
                        data = json.load(sfile)
                        meta = data.get("snapshot_meta", {})
                        snap_id = meta.get("id") or f.replace("snapshot_", "").replace(".json", "")
                        snap_type = meta.get("type", "manual" if "manual" in str(snap_id) else "eon_bounce")
                        
                        if snap_type == "eon_bounce" or str(snap_id).startswith("eon_"):
                            eon_val = meta.get("eon")
                            if eon_val is None:
                                match_e = re.search(r'eon_(\d+)', str(snap_id))
                                eon_val = int(match_e.group(1)) if match_e else "?"

                            # Strictly 1 entry per eon to prevent visual clutter in UI combobox
                            if eon_val in seen_eons:
                                continue
                            seen_eons.add(eon_val)
                            
                            sf_val = meta.get("scale_factor")
                            if sf_val is not None:
                                try:
                                    sf_f = float(sf_val)
                                    if sf_f < 6.9:
                                        norm_trans = "Rebote Gravitatorio (Túnel Cuántico)"
                                    else:
                                        norm_trans = "CCC (Muerte Térmica)"
                                except (ValueError, TypeError):
                                    raw_trans = meta.get("transition") or meta.get("label", "")
                                    norm_trans = self.normalize_transition_name(raw_trans)
                            else:
                                raw_trans = meta.get("transition") or meta.get("label", "")
                                norm_trans = self.normalize_transition_name(raw_trans)
                            
                            snapshots.append({
                                "id": str(snap_id),
                                "label": f"📷 Fin Eón {eon_val} [{norm_trans}]",
                                "type": "eon_bounce",
                                "transition": norm_trans,
                                "eon": eon_val,
                                "scale_factor": sf_val,
                                "timestamp": meta.get("timestamp", "")
                            })
                        else:
                            if snap_id in seen_manual_ids:
                                continue
                            seen_manual_ids.add(snap_id)
                            label = meta.get("label") or f"💾 Guardado: {snap_id}"
                            snapshots.append({
                                "id": str(snap_id),
                                "label": label,
                                "type": "manual",
                                "eon": meta.get("eon"),
                                "scale_factor": meta.get("scale_factor"),
                                "timestamp": meta.get("timestamp", "")
                            })
                except Exception:
                    pass

        def natural_sort_key(item):
            eon_val = item.get("eon")
            if eon_val is not None:
                try:
                    return (0, int(eon_val), item.get("timestamp", ""))
                except (ValueError, TypeError):
                    pass
            
            snap_id = str(item.get("id", ""))
            match = re.search(r'eon_(\d+)', snap_id)
            if match:
                return (0, int(match.group(1)), item.get("timestamp", ""))
            
            match_man = re.search(r'manual_(\d+)', snap_id)
            if match_man:
                return (1, int(match_man.group(1)), item.get("timestamp", ""))
            
            return (2, 0, snap_id)

        snapshots.sort(key=natural_sort_key)
        return snapshots

    def get_history_csv(self):
        """Generates a CSV string containing all historical eon metrics."""
        history = self.get_history()
        output = io.StringIO()
        writer = csv.writer(output)
        
        # CSV Headers
        writer.writerow([
            "Eon",
            "TransitionMechanism",
            "FinalScaleFactor_a",
            "PeakEntropy_kB",
            "BekensteinLimit_kB",
            "CoreMassFraction_pct",
            "FossilOdometerTotal_s",
            "PredictedH0_kms_mpc",
            "CMB_C2_Quadrupole",
            "CMB_C3_Octupole",
            "CMB_C2_C3_Ratio",
            "Planck_Chi2",
            "Attractor_X",
            "Attractor_Y",
            "Attractor_Z",
            "CPUSteps",
            "Duration_seconds",
            "Timestamp"
        ])
        
        for item in history:
            attr = item.get("attractor", {})
            obs = item.get("observational", {})
            writer.writerow([
                item.get("eon"),
                item.get("transition", "Rebote Cuántico"),
                item.get("final_scale_factor"),
                item.get("peak_s_bh"),
                item.get("s_crit"),
                item.get("core_mass_fraction"),
                item.get("fossil_odometer_total"),
                item.get("h0_predicted") if item.get("h0_predicted") is not None else obs.get("h0_predicted", self.units.H0_PLANCK_BASELINE),
                item.get("quadrupole_C2") if item.get("quadrupole_C2") is not None else obs.get("quadrupole_C2", ""),
                item.get("octopole_C3") if item.get("octopole_C3") is not None else obs.get("octopole_C3", ""),
                item.get("ratio_C2_C3") if item.get("ratio_C2_C3") is not None else obs.get("ratio_C2_C3", ""),
                item.get("planck_chi2") if item.get("planck_chi2") is not None else obs.get("planck_chi2", ""),
                attr.get("x", ""),
                attr.get("y", ""),
                attr.get("z", ""),
                item.get("eon_steps"),
                item.get("walltime_seconds"),
                item.get("timestamp")
            ])
        
        return output.getvalue()

    def save_checkpoint(self, filepath=None):
        """Saves current state to compressed .npz archive."""
        g_tag = f"_g{self.grid_size}"
        if filepath is None:
            target_filepath = os.path.join(self.checkpoint_dir, f"latest{g_tag}.npz")
        else:
            target_filepath = filepath

        # Compute Poisson Gravitational Potential Phi for complete halo & metric diagnostics
        # Includes Holographic Gravity Boost for observational consistency
        delta_rho = self.rho - self.xp.mean(self.rho)
        delta_rho_eff = delta_rho + 1.0 * self.tau_eon_start
        phi_fft = -4.0 * np.pi * self.G_CONST * self.xp.fft.fftn(delta_rho_eff) / self.k2
        phi_fft[0, 0, 0] = 0.0
        phi = self.to_cpu(self.xp.real(self.xp.fft.ifftn(phi_fft))).astype(np.float32)

        data = {
            "eon": self.eon,
            "grid_size": self.grid_size,
            "scale_factor": self.scale_factor,
            "total_steps": self.total_steps,
            "t_coord": self.t_coord,
            "seed": self.seed,
            "rng_state_json": json.dumps(self.rng.bit_generator.state),
            "rho": self.to_cpu(self.rho).astype(np.float32),
            "phi": phi,
            "T": self.to_cpu(self.T).astype(np.float32),
            "I": self.to_cpu(self.I).astype(np.float32),
            "tau": self.to_cpu(self.tau).astype(np.float32),
            "tau_physical": self.tau_physical,
            "tau_eon_start": self.to_cpu(self.tau_eon_start).astype(np.float32),
            "v_x": self.to_cpu(self.v_x).astype(np.float32),
            "v_y": self.to_cpu(self.v_y).astype(np.float32),
            "v_z": self.to_cpu(self.v_z).astype(np.float32),
            "d_tau_dt": self.to_cpu(self.d_tau_dt).astype(np.float32),
            "box_size_mpc": float(self.units.box_size_mpc),
            "h0_kms_mpc": float(self.units.h0_kms_mpc)
        }
        np.savez_compressed(target_filepath, **data)
        if filepath is None and self.grid_size == 32:
            np.savez_compressed(os.path.join(self.checkpoint_dir, "latest.npz"), **data)

    def load_checkpoint(self, filepath):
        """Loads and restores simulation state from an .npz archive."""
        if not os.path.exists(filepath):
            return False

        data = np.load(filepath)
        if "rho" in data and data["rho"].shape[0] != self.grid_size:
            return False
        self.eon = int(data["eon"])
        self.scale_factor = float(data["scale_factor"])
        self.total_steps = int(data["total_steps"])
        self.t_coord = float(data["t_coord"]) if "t_coord" in data else float(self.total_steps * self.DT)
        if "seed" in data:
            self.seed = int(data["seed"])
        if "rng_state_json" in data:
            self.rng = np.random.default_rng()
            self.rng.bit_generator.state = json.loads(str(data["rng_state_json"]))
        else:
            self.rng = np.random.default_rng(getattr(self, "seed", 42))
        self.rho = self.xp.asarray(data["rho"].astype(np.float32))
        self.T = self.xp.asarray(data["T"].astype(np.float32))
        self.I = self.xp.asarray(data["I"].astype(np.float32))
        self.tau = self.xp.asarray(data["tau"].astype(np.float32))
        self.tau_eon_start = self.xp.asarray(data["tau_eon_start"].astype(np.float32))
        self.v_x = self.xp.asarray(data["v_x"].astype(np.float32))
        self.v_y = self.xp.asarray(data["v_y"].astype(np.float32))
        self.v_z = self.xp.asarray(data["v_z"].astype(np.float32))
        self.d_tau_dt = self.xp.asarray(data["d_tau_dt"].astype(np.float32))
        self.last_bounce_step = self.total_steps
        self.eon_start_walltime = time.time()
        return True

    def _append_history(self, entry):
        """Appends a completed eon summary to the history.json log."""
        history = self.get_history()
        history.append(entry)
        with open(self.history_file, 'w', encoding='utf-8') as f:
            json.dump(history, f, indent=2)

    def get_history(self):
        """Reads and returns the list of historical eon records with normalized transition names."""
        if os.path.exists(self.history_file):
            try:
                with open(self.history_file, 'r', encoding='utf-8') as f:
                    history = json.load(f)
                    for item in history:
                        item["transition"] = self.normalize_transition_name(item.get("transition"))
                    return history
            except Exception:
                return []
        return []
