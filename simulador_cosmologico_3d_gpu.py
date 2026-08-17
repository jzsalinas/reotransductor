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

# =====================================================================
# SELECCIÓN AUTOMÁTICA DE BACKEND: GPU (CuPy / CUDA) O CPU (NumPy)
# =====================================================================
try:
    import cupy as xp
    GPU_AVAILABLE = True
    device_name = xp.cuda.runtime.getDeviceProperties(0)['name'].decode()
    BACKEND_DESC = f"GPU ({device_name})"
except Exception:
    import numpy as xp
    GPU_AVAILABLE = False
    BACKEND_DESC = "CPU (NumPy Fallback)"

def to_cpu(arr):
    """Convierte arrays de GPU a NumPy para renderizado en Matplotlib."""
    return arr.get() if hasattr(arr, 'get') else arr

# =====================================================================
# CONFIGURACIÓN DEL MODELO COSMOLÓGICO 3D Y CONSTANTES ASTROFÍSICAS
# =====================================================================
GRID_SIZE = 32              # Resolución de la malla espacial 3D (32 x 32 x 32 = 32,768 celdas)
DT = 0.05                   # Paso temporal de integración
DIFFUSION_COEFF = 0.3       # Conductividad térmica del plasma intergaláctico (k)
KAPPA = 50.0                # Constante de acoplamiento del Reotransductor (κ)
LANDAUER_DECAY = 0.015      # Tasa de decaimiento entrópico de Landauer (γ)
G_CONST = 0.04              # Constante de gravitación efectiva 3D normalizada (G)
H_0 = 0.0003                # Tasa de expansión asintótica de Hubble
CS2 = 0.18                  # Velocidad del sonido al cuadrado (Presión de Jeans / Estabilización)
C_LIGHT = 2.5               # Velocidad de la luz y límite causal relativista (v <= c)
INFLATION_BOOST = 8.0       # Impulso de super-expansión inflacionaria primordial

# Termodinámica Cósmica Real (Escala Kelvin Astrofísica):
T_CMB = 2.7255              # Temperatura de fondo cósmico real (Fondo CMB medido por COBE/Planck)
T_RECOMB = 3000.0           # Temperatura de ionización y recombinación del hidrógeno primordial (~3000 K)
GAMMA_ADIABATIC = 5.0 / 3.0 # Índice adiabático para plasma monoatómico ideal (γ = 5/3, γ-1 = 2/3)

# Capacidad holográfica volumétrica 3D calibrada a M0 = 120.0
M0_3D_CORE = 120.0
ZETA_BEKENSTEIN = 3500.0

# 1. Espectro de perturbaciones primordiales 3D P(k) ~ k^(-0.75)
kx_np = 2 * np.pi * np.fft.fftfreq(GRID_SIZE)[:, None, None].astype(np.float32)
ky_np = 2 * np.pi * np.fft.fftfreq(GRID_SIZE)[None, :, None].astype(np.float32)
kz_np = 2 * np.pi * np.fft.fftfreq(GRID_SIZE)[None, None, :].astype(np.float32)
k2_np = kx_np**2 + ky_np**2 + kz_np**2
k2_np[0, 0, 0] = 1.0  # Evitar división por cero en modo k=0
p_k_np = 1.0 / (k2_np**0.75)
p_k_np[0, 0, 0] = 0.0

sigma_g = 2.2
gaussian_k_3d_np = np.exp(-0.5 * k2_np * (sigma_g**2)).astype(np.float32)

# Transferencia de constantes espectrales a GPU/CPU
k2 = xp.asarray(k2_np, dtype=xp.float32)
p_k = xp.asarray(p_k_np, dtype=xp.float32)
gaussian_k_3d = xp.asarray(gaussian_k_3d_np, dtype=xp.float32)

# =====================================================================
# GARANTÍA DE IDENTIDAD EXACTA DE DATOS INICIALES ENTRE CPU Y GPU
# =====================================================================
np.random.seed(42)
X_np, Y_np, Z_np = np.meshgrid(np.arange(GRID_SIZE), np.arange(GRID_SIZE), np.arange(GRID_SIZE), indexing='ij')

noise_fft_cpu = np.fft.fftn(np.random.randn(GRID_SIZE, GRID_SIZE, GRID_SIZE).astype(np.float32))
fluct_cpu = np.real(np.fft.ifftn(noise_fft_cpu * p_k_np))
fluct_cpu = (fluct_cpu - np.mean(fluct_cpu)) / np.std(fluct_cpu) * 0.35

seed_A_np = 2.8 * np.exp(-((X_np - 16.0)**2 + (Y_np - 16.0)**2 + (Z_np - 16.0)**2) / 18.0)
seed_B_np = 1.9 * np.exp(-((X_np - 24.0)**2 + (Y_np - 8.0)**2 + (Z_np - 20.0)**2) / 12.0)
void_C_np = -0.6 * np.exp(-((X_np - 8.0)**2 + (Y_np - 24.0)**2 + (Z_np - 8.0)**2) / 22.0)

rho_init_np = np.maximum(0.05, 1.0 + fluct_cpu + seed_A_np + seed_B_np + void_C_np).astype(np.float32)
T_init_np = (12.0 * (rho_init_np**0.5) + 2.73).astype(np.float32)
I_init_np = np.clip((rho_init_np - 0.5) / 2.5, 0.0, 1.0).astype(np.float32)

# Carga en memoria VRAM de la GPU
X = xp.asarray(X_np, dtype=xp.float32)
Y = xp.asarray(Y_np, dtype=xp.float32)
Z = xp.asarray(Z_np, dtype=xp.float32)

rho = xp.asarray(rho_init_np, dtype=xp.float32)
T = xp.asarray(T_init_np, dtype=xp.float32)
I = xp.asarray(I_init_np, dtype=xp.float32)
tau = xp.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=xp.float32)
tau_eon_start = xp.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=xp.float32)

v_x = xp.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=xp.float32)
v_y = xp.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=xp.float32)
v_z = xp.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=xp.float32)

scale_factor = 1.0
eon = 1
steps_per_frame = 20

# Rejilla de Alta Definición para la Esfera Celeste del CMB (Mollweide S^2)
n_lat, n_lon = 90, 180
lats = np.linspace(-np.pi / 2, np.pi / 2, n_lat)
lons = np.linspace(-np.pi, np.pi, n_lon)
LON, LAT = np.meshgrid(lons, lats)

# Vectores normales de línea de visión en VRAM
n_los_x_np = np.cos(LAT) * np.cos(LON)
n_los_y_np = np.cos(LAT) * np.sin(LON)
n_los_z_np = np.sin(LAT)

n_los_x = xp.asarray(n_los_x_np, dtype=xp.float32)
n_los_y = xp.asarray(n_los_y_np, dtype=xp.float32)
n_los_z = xp.asarray(n_los_z_np, dtype=xp.float32)

r_obs = GRID_SIZE / 2.2
cx_obs, cy_obs, cz_obs = GRID_SIZE / 2.0, GRID_SIZE / 2.0, GRID_SIZE / 2.0
coords_cmb_x = (cx_obs + r_obs * n_los_x) % GRID_SIZE
coords_cmb_y = (cy_obs + r_obs * n_los_y) % GRID_SIZE
coords_cmb_z = (cz_obs + r_obs * n_los_z) % GRID_SIZE

def sample_sphere_trilinear_gpu(arr, cx, cy, cz):
    """Muestreo trilineal continuo en la esfera en memoria VRAM con condiciones de contorno periódicas."""
    x0 = xp.floor(cx).astype(xp.int32) % GRID_SIZE
    x1 = (x0 + 1) % GRID_SIZE
    y0 = xp.floor(cy).astype(xp.int32) % GRID_SIZE
    y1 = (y0 + 1) % GRID_SIZE
    z0 = xp.floor(cz).astype(xp.int32) % GRID_SIZE
    z1 = (z0 + 1) % GRID_SIZE
    
    xd = (cx - xp.floor(cx)).astype(xp.float32)
    yd = (cy - xp.floor(cy)).astype(xp.float32)
    zd = (cz - xp.floor(cz)).astype(xp.float32)
    
    c00 = arr[x0, y0, z0] * (1.0 - xd) + arr[x1, y0, z0] * xd
    c01 = arr[x0, y0, z1] * (1.0 - xd) + arr[x1, y0, z1] * xd
    c10 = arr[x0, y1, z0] * (1.0 - xd) + arr[x1, y1, z0] * xd
    c11 = arr[x0, y1, z1] * (1.0 - xd) + arr[x1, y1, z1] * xd
    
    c0 = c00 * (1.0 - yd) + c10 * yd
    c1 = c01 * (1.0 - yd) + c11 * yd
    
    return c0 * (1.0 - zd) + c1 * zd

# =====================================================================
# OPERADORES DIFERENCIALES 3D VECTORIZADOS EN GPU
# =====================================================================
def grad_3d(f):
    """Diferencias centrales 3D periódicas ultra-rápidas."""
    gx = (xp.roll(f, -1, axis=0) - xp.roll(f, 1, axis=0)) * 0.5
    gy = (xp.roll(f, -1, axis=1) - xp.roll(f, 1, axis=1)) * 0.5
    gz = (xp.roll(f, -1, axis=2) - xp.roll(f, 1, axis=2)) * 0.5
    return gx, gy, gz

def laplacian_3d(f):
    """Laplaciano discreto 3D de 7 puntos."""
    return (
        xp.roll(f, 1, axis=0) + xp.roll(f, -1, axis=0) +
        xp.roll(f, 1, axis=1) + xp.roll(f, -1, axis=1) +
        xp.roll(f, 1, axis=2) + xp.roll(f, -1, axis=2) - 6.0 * f
    )

def trigger_white_hole_eon_3d(rho_current, tau_current):
    """
    Transición Agujero Negro -> Agujero Blanco 3D (Rovelli / LQC).
    Expulsa la masa comprimida en una onda de choque esférica tridimensional con
    inercia cinética de Hubble, recalentamiento térmico primordial y acople a la red fósil tau.
    """
    max_idx = int(to_cpu(xp.argmax(rho_current)))
    x0, y0, z0 = np.unravel_index(max_idx, (GRID_SIZE, GRID_SIZE, GRID_SIZE))
    
    dx = (X - x0 + GRID_SIZE / 2) % GRID_SIZE - GRID_SIZE / 2
    dy = (Y - y0 + GRID_SIZE / 2) % GRID_SIZE - GRID_SIZE / 2
    dz = (Z - z0 + GRID_SIZE / 2) % GRID_SIZE - GRID_SIZE / 2
    r = xp.sqrt(dx**2 + dy**2 + dz**2)
    r_safe = xp.maximum(0.8, r)
    
    # 1. Acople fósil isotrópico 3D vía Fourier
    tau_fft = xp.fft.fftn(tau_current)
    tau_smooth = xp.real(xp.fft.ifftn(tau_fft * gaussian_k_3d))
    tau_norm = (tau_smooth - xp.mean(tau_smooth)) / max(1e-4, float(to_cpu(xp.std(tau_smooth))))
    fossil_fluct = xp.clip(tau_norm * 0.4, -0.8, 0.8)
    
    # 2. Nuevas perturbaciones cuánticas primordiales 3D
    noise_fft_new = xp.fft.fftn(xp.asarray(np.random.randn(GRID_SIZE, GRID_SIZE, GRID_SIZE).astype(np.float32)))
    fluct_new = xp.real(xp.fft.ifftn(noise_fft_new * p_k))
    fluct_new = (fluct_new - xp.mean(fluct_new)) / float(to_cpu(xp.std(fluct_new))) * 0.35
    
    # 3. Plasma expulsado del Agujero Blanco 3D
    core_blast = 1.8 * xp.exp(-(r**2) / 30.0)
    shell_blast = 1.5 * xp.exp(-((r - 9.0)**2) / 18.0)
    
    rho_new = xp.maximum(0.05, 1.0 + fluct_new + fossil_fluct + core_blast + shell_blast)
    
    # 4. Velocidades de eyección sublumínica 3D (v <= c)
    v_blast = 1.8 * xp.exp(-((r - 6.0)**2) / 22.0)
    vx_blast = v_blast * (dx / r_safe)
    vy_blast = v_blast * (dy / r_safe)
    vz_blast = v_blast * (dz / r_safe)
    
    # 5. Recalentamiento Térmico Primordial
    T_new = 25.0 * (rho_new**0.5) + 35.0
    
    return rho_new, vx_blast, vy_blast, vz_blast, T_new

# =====================================================================
# INTEGRACIÓN DINÁMICA POR LOTES EN GPU (BATCH COMPUTATION)
# =====================================================================
def update_cosmology_3d_batch(rho, T, I, tau, v_x, v_y, v_z, scale_factor, steps):
    """
    Ejecuta 'steps' pasos de simulación completamente en VRAM de la GPU
    sin transferencias intermedias hacia el host CPU.
    """
    for _ in range(steps):
        # a. Expansión cósmica e Inflación 3D
        if scale_factor < 1.05:
            H_eff = H_0 * (1.0 + INFLATION_BOOST * np.exp(-(scale_factor - 1.0) / 0.015))
        else:
            H_eff = H_0
        scale_factor += H_eff * DT
        
        # b. Gravedad 3D: Ecuación de Poisson ∇²Φ = 4πG(ρ - ρ̄) vía FFT 3D
        delta_rho = rho - xp.mean(rho)
        phi_fft = -4.0 * np.pi * G_CONST * xp.fft.fftn(delta_rho) / k2
        phi_fft[0, 0, 0] = 0.0
        phi = xp.real(xp.fft.ifftn(phi_fft))
        
        grad_phi_x, grad_phi_y, grad_phi_z = grad_3d(phi)
        
        # c. Gradiente de Presión de Jeans 3D (Polítropo P = CS2 * ρ^1.3)
        P = CS2 * (rho**1.3)
        grad_P_x, grad_P_y, grad_P_z = grad_3d(P)
        
        # Aceleración neta 3D
        acc_x = -grad_phi_x - (grad_P_x / (rho + 0.2))
        acc_y = -grad_phi_y - (grad_P_y / (rho + 0.2))
        acc_z = -grad_phi_z - (grad_P_z / (rho + 0.2))
        
        # Evolución hidrodinámica 3D con inercia cinética y amortiguamiento suave de Hubble
        v_x = 0.92 * v_x + 0.08 * (0.06 * acc_x)
        v_y = 0.92 * v_y + 0.08 * (0.06 * acc_y)
        v_z = 0.92 * v_z + 0.08 * (0.06 * acc_z)
        
        # Límite Causal Relativista (v <= c)
        v_mag = xp.sqrt(v_x**2 + v_y**2 + v_z**2)
        v_limit = xp.maximum(1.0, v_mag / C_LIGHT)
        v_x /= v_limit
        v_y /= v_limit
        v_z /= v_limit
        
        # d. Hidrodinámica de materia cósmica 3D
        div_flux = grad_3d(rho * v_x)[0] + grad_3d(rho * v_y)[1] + grad_3d(rho * v_z)[2]
        lap_rho = laplacian_3d(rho)
        rho = xp.clip(rho - div_flux * DT + 0.04 * lap_rho * DT, 0.02, 12.0)
        
        # e. Evolución térmica 3D: Conducción + Compresión de Choque - Enfriamiento de Hubble
        lap_T = laplacian_3d(T)
        compression_heating = 6.0 * xp.maximum(0.0, -div_flux)
        hubble_cooling = H_eff * T
        T = xp.clip(T + (DIFFUSION_COEFF * lap_T + compression_heating - hubble_cooling) * DT, 2.73, 2000.0)
        
        # f. Producción de Entropía de Onsager 3D
        inv_T = 1.0 / T
        grad_inv_T_x, grad_inv_T_y, grad_inv_T_z = grad_3d(inv_T)
        grad_T_x, grad_T_y, grad_T_z = grad_3d(T)
        
        J_T_x = -DIFFUSION_COEFF * grad_T_x
        J_T_y = -DIFFUSION_COEFF * grad_T_y
        J_T_z = -DIFFUSION_COEFF * grad_T_z
        
        sigma_thermal = xp.maximum(0.0, J_T_x * grad_inv_T_x + J_T_y * grad_inv_T_y + J_T_z * grad_inv_T_z)
        sigma_grav = (rho * (grad_phi_x**2 + grad_phi_y**2 + grad_phi_z**2)) / (T * 50.0)
        sigma_total = sigma_thermal + sigma_grav
        
        # g. Reotransductor: Tiempo Propio Emergente 3D
        d_tau_dt = KAPPA * sigma_total
        tau += d_tau_dt * DT
        
        # h. Campo Informacional Continuo 3D I(r, t)
        div_flux_I = grad_3d(I * v_x)[0] + grad_3d(I * v_y)[1] + grad_3d(I * v_z)[2]
        lap_I = laplacian_3d(I)
        sustenance = 0.6 * sigma_total * (rho / xp.mean(rho))
        thermal_noise = 0.0004 * T
        dI_dt = -div_flux_I + 0.02 * lap_I + (sustenance - thermal_noise - LANDAUER_DECAY * I)
        I = xp.clip(I + dI_dt * DT, 0.0, 1.0)
        
    return rho, T, I, tau, v_x, v_y, v_z, scale_factor, d_tau_dt, phi

def evaluate_bekenstein_trigger_3d(rho_current, tau_current, tau_eon_start, v_x, v_y, v_z, scale_factor, eon, phi_current):
    """Evalúa la condición física de rebote de Bekenstein-Hawking en la red 3D."""
    tau_current_eon = tau_current - tau_eon_start
    
    total_mass = float(to_cpu(xp.sum(rho_current)))
    core_mask = rho_current > 1.0
    core_mass = float(to_cpu(xp.sum(rho_current[core_mask])))
    mass_fraction = core_mass / max(1.0, total_mass)
    
    s_bh_eon = float(to_cpu(xp.max(tau_current_eon)))
    tau_bekenstein_crit = ZETA_BEKENSTEIN * max(1.0, (core_mass / 5000.0)**2)
    
    # Progreso unificado: condensación de masa y saturación entrópica sincronizados
    p_mass = min(1.0, mass_fraction / 0.18)
    p_entropy = min(1.0, s_bh_eon / max(1.0, tau_bekenstein_crit))
    tunnel_progress = min(1.0, min(p_mass, p_entropy))
    
    # En 3D, una concentración del 18%+ de la masa total representa un colapso maduro del núcleo
    if mass_fraction >= 0.18 and s_bh_eon >= tau_bekenstein_crit:
        eon += 1
        scale_factor = 1.0
        tau_eon_start = tau_current.copy()
        rho_current, v_x, v_y, v_z, T_current = trigger_white_hole_eon_3d(rho_current, tau_current)
        I_current = xp.clip((rho_current - 0.5) / 2.5, 0.0, 1.0)
        return rho_current, T_current, I_current, tau_current, tau_eon_start, v_x, v_y, v_z, scale_factor, eon, mass_fraction, s_bh_eon, tau_bekenstein_crit, 0.0
    
    return rho_current, None, None, tau_current, tau_eon_start, v_x, v_y, v_z, scale_factor, eon, mass_fraction, s_bh_eon, tau_bekenstein_crit, tunnel_progress

# =====================================================================
# =====================================================================
# DASHBOARD CIENTÍFICO 3D (3x3 - 9 PANELES UNIFICADOS):
# RED CÓSMICA 3D, CMB MOLLWEIDE, TELEMETRÍA, TÉRMICA, INFORMACIÓN Y TIEMPO
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

# Inicialización Arrays
z0 = GRID_SIZE // 2
rho_np = to_cpu(rho)
T_np = to_cpu(T)
I_np = to_cpu(I)
tau_np = to_cpu(tau)

# Panel 1: Red Cósmica 3D
ax_3d.view_init(elev=24, azim=-50)
ax_3d.dist = 7.8
ax_3d.set_box_aspect([1, 1, 1])
xs, ys, zs = np.where(rho_np > 1.4)
sc_3d = ax_3d.scatter(xs, ys, zs, c=rho_np[xs, ys, zs], cmap='magma', s=12, alpha=0.6, vmin=0.0, vmax=4.0)
ax_3d.set_title("1. Red Cósmica 3D (Masa / Filamentos)", fontweight='bold', fontsize=8.8, pad=8)
ax_3d.set_xlim(0, GRID_SIZE)
ax_3d.set_ylim(0, GRID_SIZE)
ax_3d.set_zlim(0, GRID_SIZE)
ax_3d.tick_params(pad=0.5, labelsize=7)
ax_3d.set_xlabel('X', fontsize=7.2)
ax_3d.set_ylabel('Y', fontsize=7.2)
ax_3d.set_zlabel('Z', fontsize=7.2)

# Panel 2: CMB Mollweide HD (Planck Style ΔT/T̄ con Muestreo Trilineal en VRAM y Efecto Doppler)
tau_s_init = sample_sphere_trilinear_gpu(tau, coords_cmb_x, coords_cmb_y, coords_cmb_z)
vx_s_init = sample_sphere_trilinear_gpu(v_x, coords_cmb_x, coords_cmb_y, coords_cmb_z)
vy_s_init = sample_sphere_trilinear_gpu(v_y, coords_cmb_x, coords_cmb_y, coords_cmb_z)
vz_s_init = sample_sphere_trilinear_gpu(v_z, coords_cmb_x, coords_cmb_y, coords_cmb_z)
v_los_init = (vx_s_init * n_los_x + vy_s_init * n_los_y + vz_s_init * n_los_z) / C_LIGHT
cmb_raw_init = xp.log10(1.0 + xp.maximum(0.0, tau_s_init)) + 0.4 * v_los_init
cmb_std_init = float(xp.maximum(1e-4, xp.std(cmb_raw_init)))
cmb_init = to_cpu((cmb_raw_init - float(xp.mean(cmb_raw_init))) / cmb_std_init)
im_cmb = ax_cmb.pcolormesh(LON, LAT, cmb_init, cmap=PLANCK_CMAP, shading='gouraud', vmin=-2.5, vmax=2.5)
title_cmb = ax_cmb.set_title("2. Fondo Cósmico CMB (Mollweide S² | ΔT/T̄)", fontweight='bold', fontsize=8.8, pad=8)
ax_cmb.grid(True, alpha=0.25)
cbar_cmb = fig.colorbar(im_cmb, ax=ax_cmb, orientation='horizontal', fraction=0.046, pad=0.10)
cbar_cmb.set_label("Anisotropías Relativas ΔT/T̄ (±2.5σ)", fontsize=7.2)

# Panel 3: Consola de Telemetría
ax_info.axis('off')
title_info = ax_info.text(0.5, 0.98, f"TELEMETRÍA 3D [{BACKEND_DESC}]", fontweight='bold', fontsize=9.5, ha='center', va='top', color='#0f172a')
telemetry_box = ax_info.text(
    0.5, 0.88, "", fontsize=7.4, fontfamily='monospace', va='top', ha='center',
    bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8fafc', edgecolor='#cbd5e1')
)

# Panel 4: Densidad de Materia Cósmica ρ
im_rho = ax_rho.imshow(rho_np[:, :, z0], cmap='magma', origin='lower', vmin=0.0, vmax=4.0)
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
im_index = ax_index.imshow(I_np[:, :, z0], cmap='cividis', origin='lower', vmin=0.0, vmax=1.0)
title_index = ax_index.set_title(f"6. Autoorganización I(r,t) (z={z0})", fontweight='bold', fontsize=8.8)
div_index = make_axes_locatable(ax_index)
cax_index = div_index.append_axes('right', size='5%', pad=0.06)
cbar_index = fig.colorbar(im_index, cax=cax_index)
cbar_index.set_label("Orden Informacional [0 a 1]", fontsize=7.2)

# Panel 7: Tiempo Emergente Lineal τ
im_tau = ax_tau.imshow(tau_np[:, :, z0], cmap='viridis', origin='lower', vmin=0.0, vmax=100.0)
title_tau = ax_tau.set_title(f"7. Tiempo Lineal (τ) (z={z0})", fontweight='bold', fontsize=8.8)
div_tau = make_axes_locatable(ax_tau)
cax_tau = div_tau.append_axes('right', size='5%', pad=0.06)
cbar_tau = fig.colorbar(im_tau, cax=cax_tau)
cbar_tau.set_label("Segundos Propios τ", fontsize=7.2)

# Panel 8: Red Fósil Logarítmica log₁₀(1+τ)
im_log_tau = ax_log_tau.imshow(np.log10(1.0 + tau_np[:, :, z0]), cmap='inferno', origin='lower', vmin=0.0, vmax=4.0)
title_log_tau = ax_log_tau.set_title(f"8. Red Fósil log₁₀(1+τ) (z={z0})", fontweight='bold', fontsize=8.8)
div_log_tau = make_axes_locatable(ax_log_tau)
cax_log_tau = div_log_tau.append_axes('right', size='5%', pad=0.06)
cbar_log_tau = fig.colorbar(im_log_tau, cax=cax_log_tau)
cbar_log_tau.set_label("Memoria Arqueológica", fontsize=7.2)

# Panel 9: Temperatura del Plasma T(r,t)
im_temp = ax_temp.imshow(T_np[:, :, z0], cmap='inferno', origin='lower', vmin=2.73, vmax=45.0)
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
    valinit=20,
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
    
    # 1. Integración física 100% en VRAM sin pausas
    rho, T, I, tau, v_x, v_y, v_z, scale_factor, d_tau_dt, phi = update_cosmology_3d_batch(
        rho, T, I, tau, v_x, v_y, v_z, scale_factor, steps_per_frame
    )
    
    # 2. Evaluación de Bekenstein al final del lote
    rho, T_new, I_new, tau, tau_eon_start, v_x, v_y, v_z, scale_factor, eon, mass_frac_val, s_bh_val, s_crit, progress = evaluate_bekenstein_trigger_3d(
        rho, tau, tau_eon_start, v_x, v_y, v_z, scale_factor, eon, phi
    )
    if T_new is not None:
        T = T_new
        I = I_new
    
    # 3. Transferencia única a CPU
    rho_cpu = to_cpu(rho)
    T_cpu = to_cpu(T)
    I_cpu = to_cpu(I)
    tau_cpu = to_cpu(tau)
    dtau_cpu = to_cpu(d_tau_dt)
    
    bx, by, bz = np.unravel_index(np.argmax(rho_cpu), rho_cpu.shape)
    z_slice = int(np.clip(bz, 0, GRID_SIZE - 1))
    
    # 4. Actualizar Scatter 3D in-place
    threshold_rho = max(1.2, float(np.percentile(rho_cpu, 88)))
    xs_n, ys_n, zs_n = np.where(rho_cpu > threshold_rho)
    if len(xs_n) > 0:
        c_vals = rho_cpu[xs_n, ys_n, zs_n]
        sc_3d._offsets3d = (xs_n, ys_n, zs_n)
        sc_3d.set_array(c_vals)
        sc_3d.set_clim(vmin=0.0, vmax=max(3.0, float(np.max(c_vals))))
    
    # 5. Actualizar Mapa CMB Mollweide HD (Memoria Fósil + Efecto Doppler Cinemático en VRAM)
    tau_s = sample_sphere_trilinear_gpu(tau, coords_cmb_x, coords_cmb_y, coords_cmb_z)
    vx_s = sample_sphere_trilinear_gpu(v_x, coords_cmb_x, coords_cmb_y, coords_cmb_z)
    vy_s = sample_sphere_trilinear_gpu(v_y, coords_cmb_x, coords_cmb_y, coords_cmb_z)
    vz_s = sample_sphere_trilinear_gpu(v_z, coords_cmb_x, coords_cmb_y, coords_cmb_z)
    v_los = (vx_s * n_los_x + vy_s * n_los_y + vz_s * n_los_z) / C_LIGHT
    cmb_raw = xp.log10(1.0 + xp.maximum(0.0, tau_s)) + 0.4 * v_los
    cmb_std = float(xp.maximum(1e-4, xp.std(cmb_raw)))
    cmb_data = to_cpu((cmb_raw - float(xp.mean(cmb_raw))) / cmb_std)
    im_cmb.set_array(cmb_data.ravel())
    im_cmb.set_clim(vmin=-2.5, vmax=2.5)
    
    # 6. Actualizar Cortes 2D en Fila 2 y Fila 3
    im_rho.set_array(rho_cpu[:, :, z_slice])
    im_rho.set_clim(vmin=0.0, vmax=max(3.0, float(np.percentile(rho_cpu[:, :, z_slice], 99))))
    
    im_rate.set_array(dtau_cpu[:, :, z_slice])
    im_rate.set_clim(vmin=0.0, vmax=max(0.5, float(np.percentile(dtau_cpu[:, :, z_slice], 99))))
    
    im_index.set_array(I_cpu[:, :, z_slice])
    
    im_tau.set_array(tau_cpu[:, :, z_slice])
    im_tau.set_clim(vmin=0.0, vmax=max(10.0, float(np.max(tau_cpu[:, :, z_slice]))))
    
    im_log_tau.set_array(np.log10(1.0 + tau_cpu[:, :, z_slice]))
    im_log_tau.set_clim(vmin=0.0, vmax=max(1.0, float(np.max(np.log10(1.0 + tau_cpu[:, :, z_slice])))))
    
    im_temp.set_array(T_cpu[:, :, z_slice])
    im_temp.set_clim(vmin=2.73, vmax=max(45.0, float(np.percentile(T_cpu[:, :, z_slice], 99))))
    
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
    fig.suptitle(f"Simulación Cosmológica 3D [GPU]: Reotransductor, CMB & Bekenstein | Eón N = {eon} [{era_str}]", fontsize=11.5, fontweight='bold')
    
    title_rho.set_text(f"4. Densidad Materia (ρ) [a={scale_factor:.2f}] (z={z_slice})")
    title_rate.set_text(f"5. Velocidad Reotransductor (dτ/dt) (z={z_slice})")
    title_index.set_text(f"6. Autoorganización I(r,t) (z={z_slice})")
    title_tau.set_text(f"7. Tiempo Lineal (τ) (z={z_slice})")
    title_log_tau.set_text(f"8. Red Fósil log₁₀(1+τ) (z={z_slice})")
    title_temp.set_text(f"9. Temperatura Plasma T (z={z_slice})")
    
    # Telemetría en consola interactiva
    bar_len = 10
    filled_len = int(progress * bar_len)
    bar_str = "█" * filled_len + "░" * (bar_len - filled_len)
    
    tau_core_val = float(np.max(tau_cpu))
    t_max_val = float(np.max(T_cpu))
    
    telemetry_str = (
        f"• Motor de Cómputo: {BACKEND_DESC}\n"
        f"• Eón Cósmico Actual: N = {eon}\n"
        f"• Factor de Escala: a = {scale_factor:.3f} (z = {redshift:.2f})\n"
        f"• Temperatura Plasma (Normalizada): {t_max_val:.1f} K ({t_max_val * 120.0:.0f} K Astrofísico)\n"
        f"• Fondo Cósmico T_CMB: 2.73 K (-270.4 °C)\n"
        f"• Velocidad de la Luz (c): {C_LIGHT:.2f} celdas/s\n"
        f"• Masa Núcleo (M_BH): {mass_frac_val * 100.0:.1f}% del total\n"
        f"• Atractor Central 3D: (x={bx}, y={by}, z={bz})\n"
        f"• Entropía Eón N (S_BH): {s_bh_val:.0f} k_B\n"
        f"• Límite Bekenstein (S_max): {s_crit:.0f} k_B\n"
        f"• Túnel Cuántico Eón {eon}: [{bar_str}] {progress * 100.0:.1f}%\n"
        f"• Odómetro Fósil Total 3D: {tau_core_val:.0f} s\n"
        f"• Estado: {'Agujero Blanco 3D Inminente' if progress >= 0.95 else 'Evolución Hidrodinámica 3D'}"
    )
    telemetry_box.set_text(telemetry_str)
    
    return sc_3d, im_cmb, im_rho, im_rate, im_index, im_tau, im_log_tau, im_temp, telemetry_box

if __name__ == '__main__':
    ani = FuncAnimation(fig, animate_3d, interval=30, blit=False, cache_frame_data=False)
    plt.show()
