import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
from matplotlib.colors import LinearSegmentedColormap

# Paleta oficial calibrada de la misión Planck (ESA Planck Legacy Archive / HEALPix)
PLANCK_CMAP = LinearSegmentedColormap.from_list(
    'planck_cmb',
    [
        '#05103a',  # Azul Marino Profundo (-2.5σ)
        '#194a8d',  # Azul Rey
        '#3288bd',  # Celeste Cian
        '#66c2a5',  # Turquesa
        '#f7f7f7',  # Neutro Blanco/Marfil (0.0σ)
        '#fee08b',  # Ámbar Claro
        '#fdae61',  # Naranja
        '#d53e4f',  # Rojo Carmesí
        '#5e001f'   # Borgoña Oscuro (+2.5σ)
    ],
    N=256
)

from server.physics_units import CosmologicalUnits, FundamentalConstants, PlanckScales

# =====================================================================
# CONFIGURACIÓN DEL MODELO COSMOLÓGICO 3D (CPU / NumPy)
# =====================================================================
GRID_SIZE = 32              # Resolución de la malla espacial 3D (32 x 32 x 32 = 32,768 celdas)
units_3d = CosmologicalUnits(box_size_mpc=100.0, grid_resolution=GRID_SIZE, c_code=2.5, h0_km_s_mpc=70.0)

DT = units_3d.DT
C_LIGHT = units_3d.c_code   # Velocidad de la luz y límite causal relativista (v <= c)
KAPPA = units_3d.get_cosmological_effective_kappa()  # Constante de acoplamiento del Reotransductor derivada
H_0 = units_3d.H_0          # Tasa de expansión asintótica de Hubble
G_CONST = units_3d.G_CONST  # Constante de gravitación efectiva 3D normalizada (G)
CS2 = units_3d.CS2_BASE     # Velocidad del sonido base a T_CMB (Presión de Jeans / Estabilización)
DIFFUSION_COEFF = units_3d.DIFFUSION_BASE  # Conductividad térmica base del plasma intergaláctico
LANDAUER_DECAY = units_3d.LANDAUER_BASE    # Tasa de decaimiento entrópico de Landauer base
INFLATION_BOOST = units_3d.INFLATION_BOOST # Impulso de super-expansión inflacionaria primordial
ZETA_BEKENSTEIN = units_3d.ZETA_BEKENSTEIN # Escala cuántica de Bekenstein-Hawking

# Termodinámica Cósmica Real (Escala Kelvin Astrofísica):
T_CMB = 2.7255              # Temperatura de fondo cósmico real (Fondo CMB medido por COBE/Planck)
T_RECOMB = 3000.0           # Temperatura de ionización y recombinación del hidrógeno primordial (~3000 K)
GAMMA_ADIABATIC = 5.0 / 3.0 # Índice adiabático para plasma monoatómico ideal (γ = 5/3, γ-1 = 2/3)

M0_3D_CORE = 120.0

# 1. Espectro de perturbaciones primordiales 3D P(k) ~ k^(-0.75)
kx = 2 * np.pi * np.fft.fftfreq(GRID_SIZE)[:, None, None].astype(np.float32)
ky = 2 * np.pi * np.fft.fftfreq(GRID_SIZE)[None, :, None].astype(np.float32)
kz = 2 * np.pi * np.fft.fftfreq(GRID_SIZE)[None, None, :].astype(np.float32)
k2 = kx**2 + ky**2 + kz**2
k2[0, 0, 0] = 1.0  # Evitar división por cero en modo k=0
p_k = 1.0 / (k2**0.75)
p_k[0, 0, 0] = 0.0

sigma_g = 2.2
gaussian_k_3d = np.exp(-0.5 * k2 * (sigma_g**2)).astype(np.float32)

# 3D Cartesian Coordinate Grids for Spatial Advection & Bounces
X, Y, Z = np.meshgrid(np.arange(GRID_SIZE), np.arange(GRID_SIZE), np.arange(GRID_SIZE), indexing='ij')

# Primordial Gaussian Random Field (GRF) with scale-invariant Harrison-Zel'dovich power spectrum P(k)
noise_fft = np.fft.fftn(np.random.randn(GRID_SIZE, GRID_SIZE, GRID_SIZE).astype(np.float32))
fluct = np.real(np.fft.ifftn(noise_fft * np.sqrt(p_k)))
fluct = (fluct - np.mean(fluct)) / max(1e-4, np.std(fluct)) * 0.45

# Organic primordial matter density field: rho_mean = 1.0 + delta_rho(k)
rho = np.maximum(0.05, 1.0 + fluct).astype(np.float32)
T = 12.0 * (rho**0.5) + 2.73
I = np.clip((rho - 0.5) / 2.5, 0.0, 1.0).astype(np.float32)
tau = np.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=np.float32)
tau_eon_start = np.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=np.float32)

v_x = np.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=np.float32)
v_y = np.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=np.float32)
v_z = np.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=np.float32)

scale_factor = 1.0
eon = 1
steps_per_frame = 5

# Rejilla de Alta Definición para la Esfera Celeste del CMB (Mollweide S^2)
n_lat, n_lon = 90, 180
lats = np.linspace(-np.pi / 2, np.pi / 2, n_lat)
lons = np.linspace(-np.pi, np.pi, n_lon)
LON, LAT = np.meshgrid(lons, lats)

# Vectores normales de línea de visión (Line-of-Sight Unit Vectors)
n_los_x = np.cos(LAT) * np.cos(LON)
n_los_y = np.cos(LAT) * np.sin(LON)
n_los_z = np.sin(LAT)

r_obs = GRID_SIZE / 2.2
cx_obs, cy_obs, cz_obs = GRID_SIZE / 2.0, GRID_SIZE / 2.0, GRID_SIZE / 2.0
coords_cmb_x = (cx_obs + r_obs * n_los_x) % GRID_SIZE
coords_cmb_y = (cy_obs + r_obs * n_los_y) % GRID_SIZE
coords_cmb_z = (cz_obs + r_obs * n_los_z) % GRID_SIZE

def sample_sphere_trilinear(arr, cx, cy, cz):
    """Muestreo trilineal continuo en la esfera con condiciones de contorno periódicas."""
    x0 = np.floor(cx).astype(int) % GRID_SIZE
    x1 = (x0 + 1) % GRID_SIZE
    y0 = np.floor(cy).astype(int) % GRID_SIZE
    y1 = (y0 + 1) % GRID_SIZE
    z0 = np.floor(cz).astype(int) % GRID_SIZE
    z1 = (z0 + 1) % GRID_SIZE
    
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

def generate_phase_locked_fluctuations(tau_current, alpha_mem=0.35):
    """Generates primordial fluctuations via Holographic Phase-Locking in Fourier Space."""
    tau_fft = np.fft.fftn(tau_current.astype(np.float32))
    theta_fossil = np.angle(tau_fft)
    
    noise_fft = np.fft.fftn(np.random.randn(GRID_SIZE, GRID_SIZE, GRID_SIZE).astype(np.float32))
    theta_quantum = np.angle(noise_fft)
    
    z_fossil = np.exp(1j * theta_fossil)
    z_quantum = np.exp(1j * theta_quantum)
    z_bounce = alpha_mem * z_fossil + (1.0 - alpha_mem) * z_quantum
    theta_bounce = np.angle(z_bounce)
    
    synthesized_fft = p_k * np.exp(1j * theta_bounce)
    synthesized_fft[0, 0, 0] = 0.0 + 0.0j
    
    fluct = np.real(np.fft.ifftn(synthesized_fft))
    std_val = float(np.std(fluct))
    if std_val > 1e-6:
        fluct = (fluct - float(np.mean(fluct))) / std_val * 0.35
    return fluct.astype(np.float32)

def trigger_white_hole_eon_3d(rho_current, tau_current):
    if np.max(rho_current) > 1.2:
        x0, y0, z0 = np.unravel_index(np.argmax(rho_current), rho_current.shape)
    else:
        x0, y0, z0 = np.unravel_index(np.argmax(tau_current), tau_current.shape)
    
    dx = (X - x0 + GRID_SIZE / 2) % GRID_SIZE - GRID_SIZE / 2
    dy = (Y - y0 + GRID_SIZE / 2) % GRID_SIZE - GRID_SIZE / 2
    dz = (Z - z0 + GRID_SIZE / 2) % GRID_SIZE - GRID_SIZE / 2
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    r_safe = np.maximum(0.8, r)
    
    # Holographic Phase-Locked Primordial Fluctuations
    fluct_new = generate_phase_locked_fluctuations(tau_current, alpha_mem=0.35)
    
    core_blast = 1.8 * np.exp(-(r**2) / 30.0)
    shell_blast = 1.5 * np.exp(-((r - 9.0)**2) / 18.0)
    
    rho_new = np.maximum(0.05, 1.0 + fluct_new + core_blast + shell_blast)
    
    v_blast = 1.8 * np.exp(-((r - 6.0)**2) / 22.0)
    vx_blast = v_blast * (dx / r_safe)
    vy_blast = v_blast * (dy / r_safe)
    vz_blast = v_blast * (dz / r_safe)
    
    # Recalentamiento Térmico Primordial
    T_new = 25.0 * (rho_new**0.5) + 35.0
    
    return rho_new, vx_blast, vy_blast, vz_blast, T_new

def update_cosmology_3d(rho, T, I, tau, tau_eon_start, v_x, v_y, v_z, scale_factor, eon):
    if scale_factor < 1.05:
        H_eff = H_0 * (1.0 + INFLATION_BOOST * np.exp(-(scale_factor - 1.0) / 0.015))
    else:
        H_eff = H_0
        
    scale_factor += H_eff * DT
    
    delta_rho = rho - np.mean(rho)
    delta_rho_fft = np.fft.fftn(delta_rho)
    phi_fft = -4.0 * np.pi * G_CONST * delta_rho_fft / k2
    phi_fft[0, 0, 0] = 0.0
    phi = np.real(np.fft.ifftn(phi_fft))
    
    grad_phi_x, grad_phi_y, grad_phi_z = np.gradient(phi)
    
    # 3. Presión de Jeans con Velocidad del Sonido Adiabática Dinámica
    cs2_field = units_3d.compute_sound_speed_sq(T, base_cs2=CS2)
    P = cs2_field * (rho**1.3)
    grad_P_x, grad_P_y, grad_P_z = np.gradient(P)
    
    acc_x = -grad_phi_x - (grad_P_x / (rho + 0.2))
    acc_y = -grad_phi_y - (grad_P_y / (rho + 0.2))
    acc_z = -grad_phi_z - (grad_P_z / (rho + 0.2))
    
    # Evolución hidrodinámica 3D con inercia cinética y amortiguamiento suave de Hubble
    v_x = 0.92 * v_x + 0.08 * (0.06 * acc_x)
    v_y = 0.92 * v_y + 0.08 * (0.06 * acc_y)
    v_z = 0.92 * v_z + 0.08 * (0.06 * acc_z)
    
    v_mag = np.sqrt(v_x**2 + v_y**2 + v_z**2)
    v_limit = np.maximum(1.0, v_mag / C_LIGHT)
    v_x /= v_limit
    v_y /= v_limit
    v_z /= v_limit
    
    flux_x = rho * v_x
    flux_y = rho * v_y
    flux_z = rho * v_z
    div_flux = np.gradient(flux_x, axis=0) + np.gradient(flux_y, axis=1) + np.gradient(flux_z, axis=2)
    
    laplacian_rho = (
        np.roll(rho, 1, axis=0) + np.roll(rho, -1, axis=0) +
        np.roll(rho, 1, axis=1) + np.roll(rho, -1, axis=1) +
        np.roll(rho, 1, axis=2) + np.roll(rho, -1, axis=2) - 6 * rho
    )
    rho_next = np.clip(rho - div_flux * DT + 0.04 * laplacian_rho * DT, 0.02, 12.0)
    
    # 5. Difusión Térmica de Spitzer-Braginskii en Plasma Cósmico
    kappa_spitzer = units_3d.compute_spitzer_conductivity(T, rho, base_k=DIFFUSION_COEFF)
    laplacian_T = (
        np.roll(T, 1, axis=0) + np.roll(T, -1, axis=0) +
        np.roll(T, 1, axis=1) + np.roll(T, -1, axis=1) +
        np.roll(T, 1, axis=2) + np.roll(T, -1, axis=2) - 6 * T
    )
    compression_heating = 6.0 * np.maximum(0.0, -div_flux)
    hubble_cooling = H_eff * T
    T_next = np.clip(T + (kappa_spitzer * laplacian_T + compression_heating - hubble_cooling) * DT, 2.73, 2000.0)
    
    inv_T = 1.0 / T_next
    grad_inv_T_x, grad_inv_T_y, grad_inv_T_z = np.gradient(inv_T)
    grad_T_x, grad_T_y, grad_T_z = np.gradient(T_next)
    
    J_T_x = -kappa_spitzer * grad_T_x
    J_T_y = -kappa_spitzer * grad_T_y
    J_T_z = -kappa_spitzer * grad_T_z
    
    sigma_thermal = np.maximum(0.0, J_T_x * grad_inv_T_x + J_T_y * grad_inv_T_y + J_T_z * grad_inv_T_z)
    sigma_grav = (rho_next * (grad_phi_x**2 + grad_phi_y**2 + grad_phi_z**2)) / (T_next * 50.0)
    sigma_total = sigma_thermal + sigma_grav
    
    d_tau_dt = KAPPA * sigma_total
    tau_next = tau + d_tau_dt * DT
    tau_current_eon = tau_next - tau_eon_start
    
    flux_I_x = I * v_x
    flux_I_y = I * v_y
    flux_I_z = I * v_z
    div_flux_I = np.gradient(flux_I_x, axis=0) + np.gradient(flux_I_y, axis=1) + np.gradient(flux_I_z, axis=2)
    
    laplacian_I = (
        np.roll(I, 1, axis=0) + np.roll(I, -1, axis=0) +
        np.roll(I, 1, axis=1) + np.roll(I, -1, axis=1) +
        np.roll(I, 1, axis=2) + np.roll(I, -1, axis=2) - 6 * I
    )
    sustenance = 0.6 * sigma_total * (rho_next / np.mean(rho_next))
    landauer_erasure = units_3d.compute_landauer_decay(T_next, base_decay=LANDAUER_DECAY)
    dI_dt = -div_flux_I + 0.02 * laplacian_I + (sustenance - landauer_erasure * I)
    I_next = np.clip(I + dI_dt * DT, 0.0, 1.0)
    
    # i. CONDICIÓN FÍSICA DE REBOTE BEKENSTEIN 3D
    total_mass = float(np.sum(rho_next))
    core_mask = rho_next > 1.0
    core_mass = float(np.sum(rho_next[core_mask]))
    mass_fraction = core_mass / max(1.0, total_mass)
    
    s_bh_eon = float(np.max(tau_current_eon))
    tau_bekenstein_crit = units_3d.compute_bekenstein_entropy_limit(
        mass_core=core_mass,
        m0_ref=5000.0,
        zeta_base=ZETA_BEKENSTEIN
    )
    
    # Progreso unificado dual: colapso gravitatorio o dilución asintótica conforme
    p_mass = min(1.0, mass_fraction / 0.18)
    p_entropy = min(1.0, s_bh_eon / max(1.0, tau_bekenstein_crit))
    p_grav = min(p_mass, p_entropy)
    p_conformal = max(0.0, min(1.0, (scale_factor - 1.0) / 6.0))
    tunnel_progress = min(1.0, max(p_grav, p_conformal))
    
    # Disparo por Colapso Gravitatorio Maduro O por Transición Conforme de Penrose (Muerte Térmica a >= 7.0)
    if (mass_fraction >= 0.18 and s_bh_eon >= tau_bekenstein_crit) or (scale_factor >= 7.0):
        eon += 1
        scale_factor = 1.0
        tau_eon_start = tau_next.copy()
        rho_next, v_x, v_y, v_z, T_next = trigger_white_hole_eon_3d(rho_next, tau_next)
        I_next = np.clip((rho_next - 0.5) / 2.5, 0.0, 1.0)
        tunnel_progress = 0.0
    
    return rho_next, T_next, I_next, tau_next, tau_eon_start, v_x, v_y, v_z, scale_factor, eon, d_tau_dt, phi, mass_fraction, s_bh_eon, tau_bekenstein_crit, tunnel_progress

# =====================================================================
# DASHBOARD CIENTÍFICO 3D (3x3 - 9 PANELES UNIFICADOS)
# =====================================================================
fig = plt.figure(figsize=(13.6, 9.8))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.38, wspace=0.32, top=0.93, bottom=0.08, left=0.06, right=0.94)

# --- Fila 1: Macrocosmos 3D, CMB y Telemetría ---
ax_3d      = fig.add_subplot(gs[0, 0], projection='3d')
ax_cmb     = fig.add_subplot(gs[0, 1], projection='mollweide')
ax_info    = fig.add_subplot(gs[0, 2])

# --- Fila 2: Densidad, Velocidad Reotransductor y Autoorganización (Corte z) ---
ax_rho     = fig.add_subplot(gs[1, 0])
ax_rate    = fig.add_subplot(gs[1, 1])
ax_index   = fig.add_subplot(gs[1, 2])

# --- Fila 3: Tiempo Lineal, Red Fósil y Temperatura (Corte z) ---
ax_tau     = fig.add_subplot(gs[2, 0])
ax_log_tau = fig.add_subplot(gs[2, 1])
ax_temp    = fig.add_subplot(gs[2, 2])

z0 = GRID_SIZE // 2

# Panel 1: Red Cósmica 3D
ax_3d.view_init(elev=24, azim=-50)
ax_3d.dist = 7.8
ax_3d.set_box_aspect([1, 1, 1])
xs, ys, zs = np.where(rho > 1.4)
sc_3d = ax_3d.scatter(xs, ys, zs, c=rho[xs, ys, zs], cmap='magma', s=12, alpha=0.6, vmin=0.0, vmax=4.0)
ax_3d.set_title("1. Red Cósmica 3D (Masa / Filamentos)", fontweight='bold', fontsize=8.8, pad=8)
ax_3d.set_xlim(0, GRID_SIZE)
ax_3d.set_ylim(0, GRID_SIZE)
ax_3d.set_zlim(0, GRID_SIZE)
ax_3d.tick_params(pad=0.5, labelsize=7)
ax_3d.set_xlabel('X', fontsize=7.2)
ax_3d.set_ylabel('Y', fontsize=7.2)
ax_3d.set_zlabel('Z', fontsize=7.2)

# Panel 2: CMB Mollweide HD (Planck Style ΔT/T̄ con Muestreo Trilineal y Efecto Doppler)
tau_s_init = sample_sphere_trilinear(tau, coords_cmb_x, coords_cmb_y, coords_cmb_z)
vx_s_init = sample_sphere_trilinear(v_x, coords_cmb_x, coords_cmb_y, coords_cmb_z)
vy_s_init = sample_sphere_trilinear(v_y, coords_cmb_x, coords_cmb_y, coords_cmb_z)
vz_s_init = sample_sphere_trilinear(v_z, coords_cmb_x, coords_cmb_y, coords_cmb_z)
v_los_init = (vx_s_init * n_los_x + vy_s_init * n_los_y + vz_s_init * n_los_z) / C_LIGHT
cmb_raw = np.log10(1.0 + np.maximum(0.0, tau_s_init)) + 0.4 * v_los_init
cmb_std = max(1e-4, float(np.std(cmb_raw)))
cmb_init = (cmb_raw - np.mean(cmb_raw)) / cmb_std
im_cmb = ax_cmb.pcolormesh(LON, LAT, cmb_init, cmap=PLANCK_CMAP, shading='gouraud', vmin=-2.5, vmax=2.5)
title_cmb = ax_cmb.set_title("2. Fondo Cósmico CMB (Mollweide S² | ΔT/T̄)", fontweight='bold', fontsize=8.8, pad=8)
ax_cmb.grid(True, alpha=0.25)
cbar_cmb = fig.colorbar(im_cmb, ax=ax_cmb, orientation='horizontal', fraction=0.046, pad=0.10)
cbar_cmb.set_label("Anisotropías Relativas ΔT/T̄ (±2.5σ)", fontsize=7.2)

# Panel 3: Consola de Telemetría
ax_info.axis('off')
title_info = ax_info.text(0.5, 0.98, "TELEMETRÍA 3D [CPU (NumPy)]", fontweight='bold', fontsize=9.5, ha='center', va='top', color='#0f172a')
telemetry_box = ax_info.text(
    0.5, 0.88, "", fontsize=7.4, fontfamily='monospace', va='top', ha='center',
    bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8fafc', edgecolor='#cbd5e1')
)

# Panel 4: Densidad de Materia Cósmica ρ
im_rho = ax_rho.imshow(rho[:, :, z0], cmap='magma', origin='lower', vmin=0.0, vmax=4.0)
title_rho = ax_rho.set_title(f"4. Densidad Materia (ρ) (z={z0})", fontweight='bold', fontsize=8.8)
div_rho = make_axes_locatable(ax_rho)
cax_rho = div_rho.append_axes('right', size='5%', pad=0.06)
cbar_rho = fig.colorbar(im_rho, cax=cax_rho)
cbar_rho.set_label("Densidad [ρ/ρ̄]", fontsize=7.2)

# Panel 5: Velocidad del Reotransductor dτ/dt
im_rate = ax_rate.imshow(np.zeros((GRID_SIZE, GRID_SIZE)), cmap='plasma', origin='lower', vmin=0.0, vmax=1.5)
title_rate = ax_rate.set_title(f"5. Velocidad Reotransductor (dτ/dt) (z={z0})", fontweight='bold', fontsize=8.8)
div_rate = make_axes_locatable(ax_rate)
cax_rate = div_rate.append_axes('right', size='5%', pad=0.06)
cbar_rate = fig.colorbar(im_rate, cax=cax_rate)
cbar_rate.set_label("Flujo Local de Presente", fontsize=7.2)

# Panel 6: Autoorganización / Negentropía I(r,t)
im_index = ax_index.imshow(I[:, :, z0], cmap='cividis', origin='lower', vmin=0.0, vmax=1.0)
title_index = ax_index.set_title(f"6. Autoorganización I(r,t) (z={z0})", fontweight='bold', fontsize=8.8)
div_index = make_axes_locatable(ax_index)
cax_index = div_index.append_axes('right', size='5%', pad=0.06)
cbar_index = fig.colorbar(im_index, cax=cax_index)
cbar_index.set_label("Orden Informacional [0 a 1]", fontsize=7.2)

# Panel 7: Tiempo Emergente Lineal τ
im_tau = ax_tau.imshow(tau[:, :, z0], cmap='viridis', origin='lower', vmin=0.0, vmax=100.0)
title_tau = ax_tau.set_title(f"7. Tiempo Lineal (τ) (z={z0})", fontweight='bold', fontsize=8.8)
div_tau = make_axes_locatable(ax_tau)
cax_tau = div_tau.append_axes('right', size='5%', pad=0.06)
cbar_tau = fig.colorbar(im_tau, cax=cax_tau)
cbar_tau.set_label("Segundos Propios τ", fontsize=7.2)

# Panel 8: Red Fósil Logarítmica log₁₀(1+τ)
im_log_tau = ax_log_tau.imshow(np.log10(1.0 + tau[:, :, z0]), cmap='inferno', origin='lower', vmin=0.0, vmax=4.0)
title_log_tau = ax_log_tau.set_title(f"8. Red Fósil log₁₀(1+τ) (z={z0})", fontweight='bold', fontsize=8.8)
div_log_tau = make_axes_locatable(ax_log_tau)
cax_log_tau = div_log_tau.append_axes('right', size='5%', pad=0.06)
cbar_log_tau = fig.colorbar(im_log_tau, cax=cax_log_tau)
cbar_log_tau.set_label("Memoria Arqueológica", fontsize=7.2)

# Panel 9: Temperatura del Plasma T(r,t)
im_temp = ax_temp.imshow(T[:, :, z0], cmap='inferno', origin='lower', vmin=2.73, vmax=45.0)
title_temp = ax_temp.set_title(f"9. Temperatura Plasma T (z={z0})", fontweight='bold', fontsize=8.8)
div_temp = make_axes_locatable(ax_temp)
cax_temp = div_temp.append_axes('right', size='5%', pad=0.06)
cbar_temp = fig.colorbar(im_temp, cax=cax_temp)
cbar_temp.set_label("Kelvins [K] (0°C = 273.15 K)", fontsize=7.2)

# Slider de Aceleración 3D (hasta 1000x)
ax_slider = fig.add_axes([0.22, 0.022, 0.56, 0.022])
slider_speed = Slider(
    ax=ax_slider,
    label='Aceleración (Pasos / Frame) ',
    valmin=1,
    valmax=1000,
    valinit=5,
    valstep=1,
    color='#2563eb'
)

def update_speed(val):
    global steps_per_frame
    steps_per_frame = int(slider_speed.val)

slider_speed.on_changed(update_speed)

# Bucle de animación continuo 3D
def animate_3d(frame):
    global rho, T, I, tau, tau_eon_start, v_x, v_y, v_z, scale_factor, eon, sc_3d
    
    s_bh_val = 0.0
    s_crit = 0.0
    progress = 0.0
    mass_frac_val = 0.0
    
    for _ in range(steps_per_frame):
        rho, T, I, tau, tau_eon_start, v_x, v_y, v_z, scale_factor, eon, d_tau_dt, phi, mass_frac_val, s_bh_val, s_crit, progress = update_cosmology_3d(
            rho, T, I, tau, tau_eon_start, v_x, v_y, v_z, scale_factor, eon
        )
    
    bx, by, bz = np.unravel_index(np.argmax(rho), rho.shape)
    z_slice = int(np.clip(bz, 0, GRID_SIZE - 1))
    
    # 1. Actualizar Scatter 3D in-place
    threshold_rho = max(1.2, float(np.percentile(rho, 88)))
    xs_n, ys_n, zs_n = np.where(rho > threshold_rho)
    if len(xs_n) > 0:
        c_vals = rho[xs_n, ys_n, zs_n]
        sc_3d._offsets3d = (xs_n, ys_n, zs_n)
        sc_3d.set_array(c_vals)
        sc_3d.set_clim(vmin=0.0, vmax=max(3.0, float(np.max(c_vals))))
    
    # 2. Actualizar Mapa CMB Mollweide HD (Sachs-Wolfe + Plasma Primordial + Memoria Fósil + Doppler)
    tau_s = sample_sphere_trilinear(tau, coords_cmb_x, coords_cmb_y, coords_cmb_z)
    T_s = sample_sphere_trilinear(T, coords_cmb_x, coords_cmb_y, coords_cmb_z)
    rho_s = sample_sphere_trilinear(rho, coords_cmb_x, coords_cmb_y, coords_cmb_z)
    vx_s = sample_sphere_trilinear(v_x, coords_cmb_x, coords_cmb_y, coords_cmb_z)
    vy_s = sample_sphere_trilinear(v_y, coords_cmb_x, coords_cmb_y, coords_cmb_z)
    vz_s = sample_sphere_trilinear(v_z, coords_cmb_x, coords_cmb_y, coords_cmb_z)
    
    v_los = (vx_s * n_los_x + vy_s * n_los_y + vz_s * n_los_z) / C_LIGHT
    mean_T = max(1e-4, float(np.mean(T)))
    mean_rho = max(1e-4, float(np.mean(rho)))
    mean_tau = max(1e-4, float(np.mean(tau)) + 1.0)

    delta_T_int = (T_s - mean_T) / mean_T
    delta_rho = (rho_s - mean_rho) / mean_rho
    delta_tau = (tau_s - float(np.mean(tau))) / mean_tau

    cmb_raw = delta_T_int + (1.0 / 3.0) * delta_rho + v_los + 0.35 * delta_tau
    cmb_std = max(1e-4, float(np.std(cmb_raw)))
    cmb_data = (cmb_raw - float(np.mean(cmb_raw))) / cmb_std
    im_cmb.set_array(cmb_data.ravel())
    im_cmb.set_clim(vmin=-2.5, vmax=2.5)
    
    # 3. Actualizar Cortes 2D en Fila 2 y Fila 3
    im_rho.set_array(rho[:, :, z_slice])
    im_rho.set_clim(vmin=0.0, vmax=max(3.0, float(np.percentile(rho[:, :, z_slice], 99))))
    
    im_rate.set_array(d_tau_dt[:, :, z_slice])
    im_rate.set_clim(vmin=0.0, vmax=max(0.5, float(np.percentile(d_tau_dt[:, :, z_slice], 99))))
    
    im_index.set_array(I[:, :, z_slice])
    
    im_tau.set_array(tau[:, :, z_slice])
    im_tau.set_clim(vmin=0.0, vmax=max(10.0, float(np.max(tau[:, :, z_slice]))))
    
    im_log_tau.set_array(np.log10(1.0 + tau[:, :, z_slice]))
    im_log_tau.set_clim(vmin=0.0, vmax=max(1.0, float(np.max(np.log10(1.0 + tau[:, :, z_slice])))))
    
    im_temp.set_array(T[:, :, z_slice])
    im_temp.set_clim(vmin=2.73, vmax=max(45.0, float(np.percentile(T[:, :, z_slice], 99))))
    
    # Identificación de la era interna del Eón
    if scale_factor < 1.05:
        era_str = "Fase de Inflación Cuántica Primordial"
    elif scale_factor < 2.5:
        era_str = "Era de Filamentos y Panqueques 3D"
    elif scale_factor < 7.0:
        era_str = "Era de Fusiones y Acreción 3D"
    elif mass_frac_val >= 0.35:
        era_str = "Era del Agujero Negro Virializado 3D"
    else:
        era_str = "Fase Asintótica Pre-Rebote 3D"
        
    redshift = max(0.0, (1.0 / scale_factor) - 1.0)
    fig.suptitle(f"Simulación Cosmológica 3D [CPU]: Reotransductor, CMB & Bekenstein | Eón N = {eon} [{era_str}]", fontsize=11.5, fontweight='bold')
    
    title_rho.set_text(f"4. Densidad Materia (ρ) [a={scale_factor:.2f}] (z={z_slice})")
    title_rate.set_text(f"5. Velocidad Reotransductor (dτ/dt) (z={z_slice})")
    title_index.set_text(f"6. Autoorganización I(r,t) (z={z_slice})")
    title_tau.set_text(f"7. Tiempo Lineal (τ) (z={z_slice})")
    title_log_tau.set_text(f"8. Red Fósil log₁₀(1+τ) (z={z_slice})")
    title_temp.set_text(f"9. Temperatura Plasma T (z={z_slice})")
    
    bar_len = 10
    filled_len = int(progress * bar_len)
    bar_str = "█" * filled_len + "░" * (bar_len - filled_len)
    
    tau_core_val = float(np.max(tau))
    t_max_val = float(np.max(T))
    
    p_conformal = max(0.0, min(1.0, (scale_factor - 1.0) / 6.0))
    p_grav = min(mass_frac_val / 0.18, s_bh_val / max(1.0, s_crit))
    prog_label_3d = f"Frontera Conforme CCC Eón {eon}" if (p_conformal > p_grav and mass_frac_val < 0.10) else f"Túnel Cuántico Eón {eon}"

    telemetry_str = (
        f"• Motor de Cómputo: CPU (NumPy)\n"
        f"• Eón Cósmico Actual: N = {eon}\n"
        f"• Factor de Escala: a = {scale_factor:.3f} (z = {redshift:.2f})\n"
        f"• Temperatura Plasma (Normalizada): {t_max_val:.1f} K ({t_max_val * 120.0:.0f} K Astrofísico)\n"
        f"• Fondo Cósmico T_CMB: 2.73 K (-270.4 °C)\n"
        f"• Velocidad de la Luz (c): {C_LIGHT:.2f} celdas/s\n"
        f"• Masa Núcleo (M_BH): {mass_frac_val * 100.0:.1f}% del total\n"
        f"• Atractor Central 3D: (x={bx}, y={by}, z={bz})\n"
        f"• Entropía Eón N (S_BH): {s_bh_val:.0f} k_B\n"
        f"• Límite Bekenstein (S_max): {s_crit:.0f} k_B\n"
        f"• {prog_label_3d}: [{bar_str}] {progress * 100.0:.1f}%\n"
        f"• Odómetro Fósil Total 3D: {tau_core_val:.0f} s\n"
        f"• Estado: {'Agujero Blanco 3D Inminente' if progress >= 0.95 else 'Evolución Hidrodinámica 3D'}"
    )
    telemetry_box.set_text(telemetry_str)
    
    return sc_3d, im_cmb, im_rho, im_rate, im_index, im_tau, im_log_tau, im_temp, telemetry_box

if __name__ == '__main__':
    ani = FuncAnimation(fig, animate_3d, interval=30, blit=False, cache_frame_data=False)
    plt.show()
