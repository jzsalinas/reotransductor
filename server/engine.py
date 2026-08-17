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
import numpy as np
from server.notifier import TelegramNotifier

class CosmologicalEngine:
    """
    Autonomous 3D Cosmological Physics Engine.
    Executes Navier-Stokes, Poisson Gravity, Onsager Emergent Time, and Bekenstein Quantum Bounce.
    """

    def __init__(self, grid_size=32, checkpoint_dir="checkpoints", auto_resume=True, force_reset=False):
        self.grid_size = grid_size
        self.checkpoint_dir = checkpoint_dir
        self.snapshots_dir = os.path.join(self.checkpoint_dir, "snapshots")
        os.makedirs(self.checkpoint_dir, exist_ok=True)
        os.makedirs(self.snapshots_dir, exist_ok=True)
        self.history_file = os.path.join(self.checkpoint_dir, "history.json")

        # Telegram Notifier (loads from local gitignored telegram_config.json)
        self.notifier = TelegramNotifier()

        # Physical Constants (Exact match with simulador_cosmologico_3d.py)
        self.DT = 0.05
        self.DIFFUSION_COEFF = 0.3
        self.KAPPA = 50.0
        self.LANDAUER_DECAY = 0.015
        self.G_CONST = 0.04
        self.H_0 = 0.0003
        self.CS2 = 0.18
        self.C_LIGHT = 2.5
        self.INFLATION_BOOST = 8.0
        self.ZETA_BEKENSTEIN = 3500.0
        self.MASS_THRESHOLD = 0.18
        self.M0_CORE = 5000.0

        # Spatial Coordinates & Fourier Mesh
        self.X, self.Y, self.Z = np.meshgrid(
            np.arange(self.grid_size),
            np.arange(self.grid_size),
            np.arange(self.grid_size),
            indexing='ij'
        )

        kx = 2.0 * np.pi * np.fft.fftfreq(self.grid_size)[:, None, None].astype(np.float32)
        ky = 2.0 * np.pi * np.fft.fftfreq(self.grid_size)[None, :, None].astype(np.float32)
        kz = 2.0 * np.pi * np.fft.fftfreq(self.grid_size)[None, None, :].astype(np.float32)
        self.k2 = kx**2 + ky**2 + kz**2
        self.k2[0, 0, 0] = 1.0  # Regularize DC mode

        self.p_k = 1.0 / (self.k2**0.75)
        self.p_k[0, 0, 0] = 0.0

        sigma_g = 2.2
        self.gaussian_k_3d = np.exp(-0.5 * self.k2 * (sigma_g**2)).astype(np.float32)

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
        self.steps_per_frame = 20
        self.total_steps = 0
        self.eon = 1
        self.scale_factor = 1.0
        self.s_bh_val = 0.0
        self.s_crit = 3500.0
        self.progress = 0.0
        self.mass_frac_val = 0.0
        self.last_bounce_step = 0
        self.eon_start_walltime = time.time()

        # Initialize or Resume
        latest_checkpoint = os.path.join(self.checkpoint_dir, "latest.npz")
        if force_reset:
            self.reset_simulation(archive_existing=True)
        elif auto_resume and os.path.exists(latest_checkpoint):
            self.load_checkpoint(latest_checkpoint)
        else:
            self._init_primordial_state()

    def _init_primordial_state(self):
        """Initializes primordial cosmological fields identically to local 3D simulator."""
        np.random.seed(42)
        noise_fft = np.fft.fftn(np.random.randn(self.grid_size, self.grid_size, self.grid_size).astype(np.float32))
        fluct = np.real(np.fft.ifftn(noise_fft * self.p_k))
        fluct = (fluct - np.mean(fluct)) / np.std(fluct) * 0.35

        seed_A = 2.8 * np.exp(-((self.X - 16.0)**2 + (self.Y - 16.0)**2 + (self.Z - 16.0)**2) / 18.0)
        seed_B = 1.9 * np.exp(-((self.X - 24.0)**2 + (self.Y - 8.0)**2 + (self.Z - 20.0)**2) / 12.0)
        void_C = -0.6 * np.exp(-((self.X - 8.0)**2 + (self.Y - 24.0)**2 + (self.Z - 8.0)**2) / 22.0)

        self.rho = np.maximum(0.05, 1.0 + fluct + seed_A + seed_B + void_C).astype(np.float32)
        self.T = (12.0 * (self.rho**0.5) + 2.73).astype(np.float32)
        self.I = np.clip((self.rho - 0.5) / 2.5, 0.0, 1.0).astype(np.float32)
        self.tau = np.zeros((self.grid_size, self.grid_size, self.grid_size), dtype=np.float32)
        self.tau_eon_start = np.zeros((self.grid_size, self.grid_size, self.grid_size), dtype=np.float32)
        self.d_tau_dt = np.zeros((self.grid_size, self.grid_size, self.grid_size), dtype=np.float32)

        self.v_x = np.zeros((self.grid_size, self.grid_size, self.grid_size), dtype=np.float32)
        self.v_y = np.zeros((self.grid_size, self.grid_size, self.grid_size), dtype=np.float32)
        self.v_z = np.zeros((self.grid_size, self.grid_size, self.grid_size), dtype=np.float32)

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
        x0 = np.floor(cx).astype(int) % self.grid_size
        x1 = (x0 + 1) % self.grid_size
        y0 = np.floor(cy).astype(int) % self.grid_size
        y1 = (y0 + 1) % self.grid_size
        z0 = np.floor(cz).astype(int) % self.grid_size
        z1 = (z0 + 1) % self.grid_size

        xd = (cx - np.floor(cx)).astype(np.float32)
        yd = (cy - np.floor(cy)).astype(np.float32)
        zd = (cz - np.floor(cz)).astype(np.float32)

        c00 = arr[x0, y0, z0] * (1.0 - xd) + arr[x1, y0, z0] * xd
        c01 = arr[x0, y0, z1] * (1.0 - xd) + arr[x1, y0, z1] * xd
        c10 = arr[x0, y1, z0] * (1.0 - xd) + arr[x1, y1, z0] * xd
        c11 = arr[x0, y1, z1] * (1.0 - xd) + arr[x1, y1, z1] * xd

        c0 = c00 * (1.0 - yd) + c10 * yd
        c1 = c01 * (1.0 - yd) + c11 * yd

        return c0 * (1.0 - zd) + c1 * zd

    def _trigger_white_hole_eon_3d(self):
        """Detonates a white hole quantum bounce and initializes the next eon."""
        x0, y0, z0 = np.unravel_index(np.argmax(self.rho), self.rho.shape)

        dx = (self.X - x0 + self.grid_size / 2.0) % self.grid_size - self.grid_size / 2.0
        dy = (self.Y - y0 + self.grid_size / 2.0) % self.grid_size - self.grid_size / 2.0
        dz = (self.Z - z0 + self.grid_size / 2.0) % self.grid_size - self.grid_size / 2.0
        r = np.sqrt(dx**2 + dy**2 + dz**2)
        r_safe = np.maximum(0.8, r)

        tau_fft = np.fft.fftn(self.tau)
        tau_smooth = np.real(np.fft.ifftn(tau_fft * self.gaussian_k_3d))
        tau_smooth = (tau_smooth - np.mean(tau_smooth)) / max(1e-5, np.std(tau_smooth))

        noise_raw = np.fft.fftn(np.random.randn(self.grid_size, self.grid_size, self.grid_size).astype(np.float32))
        fluct_new = np.real(np.fft.ifftn(noise_raw * self.p_k))
        fluct_new = (fluct_new - np.mean(fluct_new)) / max(1e-5, np.std(fluct_new))

        primordial_blast = 3.5 * np.exp(-r**2 / 16.0)
        thermal_reheating = 85.0 * np.exp(-r**2 / 20.0) + 2.73

        v_exp_mag = 2.4 * np.exp(-r**2 / 25.0) * (r / self.grid_size)
        v_x_new = v_exp_mag * (dx / r_safe)
        v_y_new = v_exp_mag * (dy / r_safe)
        v_z_new = v_exp_mag * (dz / r_safe)

        rho_new = np.clip(1.0 + 0.35 * fluct_new + 0.25 * tau_smooth + primordial_blast, 0.05, 12.0)
        T_new = np.clip(thermal_reheating + 15.0 * np.abs(fluct_new), 2.73, 2000.0)

        return rho_new, v_x_new, v_y_new, v_z_new, T_new

    def step(self):
        """Executes a single Runge-Kutta / Eulerian cosmological differential step."""
        # 1. Cosmological Scale Factor Evolution
        if self.scale_factor < 1.05:
            H_eff = self.H_0 * (1.0 + self.INFLATION_BOOST * np.exp(-(self.scale_factor - 1.0) / 0.015))
        else:
            H_eff = self.H_0

        self.scale_factor += H_eff * self.DT

        # 2. Gravitational Potential via 3D Poisson Equation (FFT)
        delta_rho = self.rho - np.mean(self.rho)
        delta_rho_fft = np.fft.fftn(delta_rho)
        phi_fft = -4.0 * np.pi * self.G_CONST * delta_rho_fft / self.k2
        phi_fft[0, 0, 0] = 0.0
        phi = np.real(np.fft.ifftn(phi_fft))

        grad_phi_x, grad_phi_y, grad_phi_z = np.gradient(phi)

        # 3. Jeans Thermal Pressure & Navier-Stokes Acceleration
        P = self.CS2 * (self.rho**1.3)
        grad_P_x, grad_P_y, grad_P_z = np.gradient(P)

        acc_x = -grad_phi_x - (grad_P_x / (self.rho + 0.2))
        acc_y = -grad_phi_y - (grad_P_y / (self.rho + 0.2))
        acc_z = -grad_phi_z - (grad_P_z / (self.rho + 0.2))

        self.v_x = 0.92 * self.v_x + 0.08 * (0.06 * acc_x)
        self.v_y = 0.92 * self.v_y + 0.08 * (0.06 * acc_y)
        self.v_z = 0.92 * self.v_z + 0.08 * (0.06 * acc_z)

        # Relativistic Causal Limit (v <= c)
        v_mag = np.sqrt(self.v_x**2 + self.v_y**2 + self.v_z**2)
        v_limit = np.maximum(1.0, v_mag / self.C_LIGHT)
        self.v_x /= v_limit
        self.v_y /= v_limit
        self.v_z /= v_limit

        # 4. Matter Continuity Equation
        flux_x = self.rho * self.v_x
        flux_y = self.rho * self.v_y
        flux_z = self.rho * self.v_z
        div_flux = np.gradient(flux_x, axis=0) + np.gradient(flux_y, axis=1) + np.gradient(flux_z, axis=2)

        laplacian_rho = (
            np.roll(self.rho, 1, axis=0) + np.roll(self.rho, -1, axis=0) +
            np.roll(self.rho, 1, axis=1) + np.roll(self.rho, -1, axis=1) +
            np.roll(self.rho, 1, axis=2) + np.roll(self.rho, -1, axis=2) - 6.0 * self.rho
        )
        self.rho = np.clip(self.rho - div_flux * self.DT + 0.04 * laplacian_rho * self.DT, 0.02, 12.0)

        # 5. Thermal Field & Gravitational Compression Heating
        laplacian_T = (
            np.roll(self.T, 1, axis=0) + np.roll(self.T, -1, axis=0) +
            np.roll(self.T, 1, axis=1) + np.roll(self.T, -1, axis=1) +
            np.roll(self.T, 1, axis=2) + np.roll(self.T, -1, axis=2) - 6.0 * self.T
        )
        compression_heating = 6.0 * np.maximum(0.0, -div_flux)
        hubble_cooling = H_eff * self.T
        self.T = np.clip(self.T + (self.DIFFUSION_COEFF * laplacian_T + compression_heating - hubble_cooling) * self.DT, 2.73, 2000.0)

        # 6. Onsager Irreversible Entropy Production & Reotransductor Time
        inv_T = 1.0 / self.T
        grad_inv_T_x, grad_inv_T_y, grad_inv_T_z = np.gradient(inv_T)
        grad_T_x, grad_T_y, grad_T_z = np.gradient(self.T)

        J_T_x = -self.DIFFUSION_COEFF * grad_T_x
        J_T_y = -self.DIFFUSION_COEFF * grad_T_y
        J_T_z = -self.DIFFUSION_COEFF * grad_T_z

        sigma_thermal = np.maximum(0.0, J_T_x * grad_inv_T_x + J_T_y * grad_inv_T_y + J_T_z * grad_inv_T_z)
        sigma_grav = (self.rho * (grad_phi_x**2 + grad_phi_y**2 + grad_phi_z**2)) / (self.T * 50.0)
        sigma_total = sigma_thermal + sigma_grav

        self.d_tau_dt = self.KAPPA * sigma_total
        self.tau += self.d_tau_dt * self.DT

        # 7. Landauer Negentropy / Informational Field
        flux_I_x = self.I * self.v_x
        flux_I_y = self.I * self.v_y
        flux_I_z = self.I * self.v_z
        div_flux_I = np.gradient(flux_I_x, axis=0) + np.gradient(flux_I_y, axis=1) + np.gradient(flux_I_z, axis=2)

        laplacian_I = (
            np.roll(self.I, 1, axis=0) + np.roll(self.I, -1, axis=0) +
            np.roll(self.I, 1, axis=1) + np.roll(self.I, -1, axis=1) +
            np.roll(self.I, 1, axis=2) + np.roll(self.I, -1, axis=2) - 6.0 * self.I
        )
        sustenance = 0.6 * sigma_total * (self.rho / np.mean(self.rho))
        thermal_noise = 0.0004 * self.T
        dI_dt = -div_flux_I + 0.02 * laplacian_I + (sustenance - thermal_noise - self.LANDAUER_DECAY * self.I)
        self.I = np.clip(self.I + dI_dt * self.DT, 0.0, 1.0)

        # 8. Bekenstein Quantum Saturation & Eon Bounce Trigger
        tau_current_eon = self.tau - self.tau_eon_start
        total_mass = float(np.sum(self.rho))
        core_mask = self.rho > 1.0
        core_mass = float(np.sum(self.rho[core_mask]))
        self.mass_frac_val = core_mass / max(1.0, total_mass)

        self.s_bh_val = float(np.max(tau_current_eon))
        self.s_crit = self.ZETA_BEKENSTEIN * max(1.0, (core_mass / self.M0_CORE)**2)

        p_mass = min(1.0, self.mass_frac_val / self.MASS_THRESHOLD)
        p_entropy = min(1.0, self.s_bh_val / max(1.0, self.s_crit))
        self.progress = min(1.0, min(p_mass, p_entropy))

        self.total_steps += 1

        if self.mass_frac_val >= self.MASS_THRESHOLD and self.s_bh_val >= self.s_crit:
            self._handle_bounce()

    def _handle_bounce(self):
        """Processes transition to next eon, logs history, archives full snapshot, and notifies Telegram."""
        eon_duration_wall = time.time() - self.eon_start_walltime
        eon_steps = self.total_steps - self.last_bounce_step

        # Save Full Visual Snapshot of the Completed Eon at Bounce
        final_snapshot = self.get_visual_payload()
        final_snapshot["snapshot_meta"] = {
            "id": f"eon_{self.eon}",
            "label": f"📷 Fin Eón {self.eon} [Rebote Cuántico]",
            "type": "eon_bounce",
            "eon": self.eon,
            "scale_factor": round(float(self.scale_factor), 3),
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        snapshot_path = os.path.join(self.snapshots_dir, f"snapshot_eon_{self.eon}.json")
        try:
            with open(snapshot_path, "w", encoding="utf-8") as f:
                json.dump(final_snapshot, f)
        except Exception:
            pass

        # Record History Entry
        history_entry = {
            "eon": self.eon,
            "final_scale_factor": round(float(self.scale_factor), 3),
            "peak_s_bh": round(float(self.s_bh_val), 1),
            "s_crit": round(float(self.s_crit), 1),
            "core_mass_fraction": round(float(self.mass_frac_val * 100.0), 2),
            "fossil_odometer_total": round(float(np.max(self.tau)), 1),
            "eon_steps": eon_steps,
            "walltime_seconds": round(eon_duration_wall, 1),
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        self._append_history(history_entry)

        # Trigger Telegram Alert if configured
        self.notifier.check_and_notify_eon(history_entry)

        # Increment Eon
        self.eon += 1
        self.scale_factor = 1.0
        self.tau_eon_start = self.tau.copy()
        self.rho, self.v_x, self.v_y, self.v_z, self.T = self._trigger_white_hole_eon_3d()
        self.I = np.clip((self.rho - 0.5) / 2.5, 0.0, 1.0)
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
        bx, by, bz = np.unravel_index(np.argmax(self.rho), self.rho.shape)
        z_slice = int(np.clip(bz, 0, self.grid_size - 1))

        if self.scale_factor < 1.05:
            era_str = "Fase de Inflación Cuántica Primordial"
        elif self.scale_factor < 2.5:
            era_str = "Era de Filamentos y Panqueques 3D"
        elif self.scale_factor < 7.0:
            era_str = "Era de Fusiones y Acreción 3D"
        elif self.mass_frac_val >= 0.35:
            era_str = "Era del Agujero Negro Virializado 3D"
        else:
            era_str = "Fase Asintótica Pre-Rebote 3D"

        redshift = max(0.0, (1.0 / self.scale_factor) - 1.0)
        temp_norm = float(np.percentile(self.T, 99))
        temp_astro = temp_norm * 120.0

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
            "fossil_odometer": round(float(np.max(self.tau)), 0),
            "attractor": {"x": int(bx), "y": int(by), "z": int(bz)},
            "z_slice": z_slice,
            "total_steps": self.total_steps,
            "is_running": self.is_running,
            "steps_per_frame": self.steps_per_frame,
            "state_status": "Agujero Blanco 3D Inminente" if self.progress >= 0.95 else "Evolución Hidrodinámica 3D"
        }

    def get_visual_payload(self):
        """Constructs data arrays for the 9-panel web dashboard."""
        bx, by, bz = np.unravel_index(np.argmax(self.rho), self.rho.shape)
        z_slice = int(np.clip(bz, 0, self.grid_size - 1))

        # 1. 3D Cosmic Web point cloud (sparse points above threshold)
        threshold_rho = max(1.2, float(np.percentile(self.rho, 88)))
        xs, ys, zs = np.where(self.rho > threshold_rho)
        points_3d = []
        if len(xs) > 0:
            step_stride = max(1, len(xs) // 300)
            for i in range(0, len(xs), step_stride):
                points_3d.append([
                    int(xs[i]), int(ys[i]), int(zs[i]),
                    round(float(self.rho[xs[i], ys[i], zs[i]]), 2)
                ])

        # 2. HD Mollweide CMB (90x180) with Trilinear Sampling & Doppler Shift
        tau_s = self._sample_sphere_trilinear(self.tau, self.coords_cmb_x, self.coords_cmb_y, self.coords_cmb_z)
        vx_s = self._sample_sphere_trilinear(self.v_x, self.coords_cmb_x, self.coords_cmb_y, self.coords_cmb_z)
        vy_s = self._sample_sphere_trilinear(self.v_y, self.coords_cmb_x, self.coords_cmb_y, self.coords_cmb_z)
        vz_s = self._sample_sphere_trilinear(self.v_z, self.coords_cmb_x, self.coords_cmb_y, self.coords_cmb_z)

        v_los = (vx_s * self.n_los_x + vy_s * self.n_los_y + vz_s * self.n_los_z) / self.C_LIGHT
        cmb_raw = np.log10(1.0 + np.maximum(0.0, tau_s)) + 0.4 * v_los
        cmb_std = max(1e-4, float(np.std(cmb_raw)))
        cmb_norm = np.clip((cmb_raw - float(np.mean(cmb_raw))) / cmb_std, -2.5, 2.5)

        # 3. 2D Cross Sections (32x32) at z_slice
        slice_rho = np.round(self.rho[:, :, z_slice], 3).tolist()
        slice_rate = np.round(self.d_tau_dt[:, :, z_slice], 3).tolist()
        slice_index = np.round(self.I[:, :, z_slice], 3).tolist()
        slice_tau = np.round(self.tau[:, :, z_slice], 1).tolist()
        slice_log_tau = np.round(np.log10(1.0 + np.maximum(0.0, self.tau[:, :, z_slice])), 3).tolist()
        slice_temp = np.round(self.T[:, :, z_slice], 2).tolist()

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
        snapshot_id = f"manual_{timestamp_id}"
        
        payload = self.get_visual_payload()
        payload["snapshot_meta"] = {
            "id": snapshot_id,
            "label": f"💾 Guardado: Eón {self.eon} (a={self.scale_factor:.3f}, {time_display})",
            "type": "manual",
            "eon": self.eon,
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
        
        # Candidate file names
        candidates = [
            f"snapshot_{str_id}.json",
            f"snapshot_eon_{str_id}.json",
            f"snapshot_manual_{str_id}.json",
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

    def get_available_snapshots(self):
        """Returns list of all available snapshots (both automatic eon bounces and manual checkpoints)."""
        if not os.path.exists(self.snapshots_dir):
            return []
        files = os.listdir(self.snapshots_dir)
        snapshots = []
        for f in files:
            if f.startswith("snapshot_") and f.endswith(".json"):
                path = os.path.join(self.snapshots_dir, f)
                try:
                    with open(path, "r", encoding="utf-8") as sfile:
                        data = json.load(sfile)
                        meta = data.get("snapshot_meta")
                        if meta:
                            snapshots.append(meta)
                        else:
                            # Fallback for older snapshots without explicit meta
                            snap_id = f.replace("snapshot_", "").replace(".json", "")
                            if snap_id.startswith("eon_"):
                                eon_num = snap_id.replace("eon_", "")
                                snapshots.append({
                                    "id": snap_id,
                                    "label": f"📷 Fin Eón {eon_num} [Rebote Cuántico]",
                                    "type": "eon_bounce"
                                })
                            else:
                                snapshots.append({
                                    "id": snap_id,
                                    "label": f"💾 Guardado: {snap_id}",
                                    "type": "manual"
                                })
                except Exception:
                    pass
        
        # Sort by timestamp / id
        snapshots.sort(key=lambda x: x.get("id", ""), reverse=False)
        return snapshots

    def get_history_csv(self):
        """Generates a CSV string containing all historical eon metrics."""
        history = self.get_history()
        output = io.StringIO()
        writer = csv.writer(output)
        
        # CSV Headers
        writer.writerow([
            "Eon",
            "FinalScaleFactor_a",
            "PeakEntropy_kB",
            "BekensteinLimit_kB",
            "CoreMassFraction_pct",
            "FossilOdometerTotal_s",
            "CPUSteps",
            "Duration_seconds",
            "Timestamp"
        ])
        
        for item in history:
            writer.writerow([
                item.get("eon"),
                item.get("final_scale_factor"),
                item.get("peak_s_bh"),
                item.get("s_crit"),
                item.get("core_mass_fraction"),
                item.get("fossil_odometer_total"),
                item.get("eon_steps"),
                item.get("walltime_seconds"),
                item.get("timestamp")
            ])
        
        return output.getvalue()

    def save_checkpoint(self, filepath=None):
        """Saves current state to compressed .npz archive."""
        if filepath is None:
            eon_filepath = os.path.join(self.checkpoint_dir, f"eon_{self.eon}.npz")
            latest_filepath = os.path.join(self.checkpoint_dir, "latest.npz")
        else:
            eon_filepath = filepath
            latest_filepath = filepath

        data = {
            "eon": self.eon,
            "scale_factor": self.scale_factor,
            "total_steps": self.total_steps,
            "rho": self.rho,
            "T": self.T,
            "I": self.I,
            "tau": self.tau,
            "tau_eon_start": self.tau_eon_start,
            "v_x": self.v_x,
            "v_y": self.v_y,
            "v_z": self.v_z,
            "d_tau_dt": self.d_tau_dt
        }
        np.savez_compressed(eon_filepath, **data)
        if filepath is None:
            np.savez_compressed(latest_filepath, **data)

    def load_checkpoint(self, filepath):
        """Loads and restores simulation state from an .npz archive."""
        if not os.path.exists(filepath):
            return False

        data = np.load(filepath)
        self.eon = int(data["eon"])
        self.scale_factor = float(data["scale_factor"])
        self.total_steps = int(data["total_steps"])
        self.rho = data["rho"].astype(np.float32)
        self.T = data["T"].astype(np.float32)
        self.I = data["I"].astype(np.float32)
        self.tau = data["tau"].astype(np.float32)
        self.tau_eon_start = data["tau_eon_start"].astype(np.float32)
        self.v_x = data["v_x"].astype(np.float32)
        self.v_y = data["v_y"].astype(np.float32)
        self.v_z = data["v_z"].astype(np.float32)
        self.d_tau_dt = data["d_tau_dt"].astype(np.float32)
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
        """Reads and returns the list of historical eon records."""
        if os.path.exists(self.history_file):
            try:
                with open(self.history_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except Exception:
                return []
        return []
