import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider

# =====================================================================
# CONFIGURACIÓN DEL MODELO COSMOLÓGICO CÍCLICO Y PARÁMETROS FÍSICOS
# =====================================================================
GRID_SIZE = 50              # Resolución de la malla espacial cosmológica (N x N)
DT = 0.05                   # Paso temporal de integración
DIFFUSION_COEFF = 0.3       # Conductividad térmica del plasma intergaláctico (k)
KAPPA = 50.0                # Constante de acoplamiento del Reotransductor (κ)
LANDAUER_DECAY = 0.015      # Tasa de decaimiento entrópico de Landauer (γ)
G_CONST = 0.04              # Constante de gravitación efectiva normalizada (G)
H_0 = 0.0003                # Tasa de expansión asintótica de Hubble
CS2 = 0.18                  # Velocidad del sonido al cuadrado (Presión de Jeans / Estabilización)
C_LIGHT = 2.5               # Velocidad de la luz y límite causal relativista (v <= c)
INFLATION_BOOST = 8.0       # Impulso de super-expansión inflacionaria primordial al inicio de cada Eón
ZETA_BEKENSTEIN = 3500.0    # Constante de capacidad entrópica de Bekenstein calibrada (S_max ~ ZETA * (M_BH/30)^2)

# 1. Espectro de perturbaciones primordiales P(k) ~ k^(-0.6)
kx = 2 * np.pi * np.fft.fftfreq(GRID_SIZE)[:, None]
ky = 2 * np.pi * np.fft.fftfreq(GRID_SIZE)[None, :]
k2 = kx**2 + ky**2
k2[0, 0] = 1.0  # Evitar división por cero en modo k=0
p_k = 1.0 / (k2**0.6)
p_k[0, 0] = 0.0

# Filtro gaussiano espectral en Fourier para suavizado isotrópico
sigma_g = 2.5
gaussian_k = np.exp(-0.5 * k2 * (sigma_g**2))

np.random.seed(42)
X, Y = np.meshgrid(np.arange(GRID_SIZE), np.arange(GRID_SIZE))

noise_fft = np.fft.fft2(np.random.randn(GRID_SIZE, GRID_SIZE))
fluct = np.real(np.fft.ifft2(noise_fft * p_k))
fluct = (fluct - np.mean(fluct)) / np.std(fluct) * 0.3

# Semillas originales para el Eón 1:
seed_A = 2.8 * np.exp(-((X - 23.5)**2 + (Y - 23.5)**2) / 14.0)   # Supercúmulo central A
seed_B = 1.9 * np.exp(-((X - 37.5)**2 + (Y - 13.5)**2) / 10.0)   # Filamento activo B
void_C = -0.6 * np.exp(-((X - 9.5)**2 + (Y - 37.5)**2) / 18.0)   # Vacío cósmico C

# Estado primordial inicial del Eón 1
rho = np.maximum(0.05, 1.0 + fluct + seed_A + seed_B + void_C)
T = 12.0 * (rho**0.5) + 2.73
I = np.clip((rho - 0.5) / 2.5, 0.0, 1.0)
tau = np.zeros((GRID_SIZE, GRID_SIZE))
tau_eon_start = np.zeros((GRID_SIZE, GRID_SIZE))

v_x = np.zeros((GRID_SIZE, GRID_SIZE))
v_y = np.zeros((GRID_SIZE, GRID_SIZE))

scale_factor = 1.0
eon = 1
steps_per_frame = 5

def trigger_white_hole_eon(rho_current, tau_current):
    """
    Transición Agujero Negro -> Agujero Blanco & Inflación Cuántica Post-Rebote (Rovelli / LQC / Guth).
    Expulsa la materia acumulada en una onda de choque radial sublumínica (v <= c),
    recalienta el plasma (reheating) y acopla las nuevas semillas cósmicas a la memoria fósil tau.
    """
    y0, x0 = np.unravel_index(np.argmax(rho_current), rho_current.shape)
    
    dx = (X - x0 + GRID_SIZE / 2) % GRID_SIZE - GRID_SIZE / 2
    dy = (Y - y0 + GRID_SIZE / 2) % GRID_SIZE - GRID_SIZE / 2
    r = np.sqrt(dx**2 + dy**2)
    r_safe = np.maximum(0.8, r)
    
    # 1. Acople fósil isotrópico suave vía Fourier (elimina artefactos en cruz)
    tau_fft = np.fft.fft2(tau_current)
    tau_smooth = np.real(np.fft.ifft2(tau_fft * gaussian_k))
    tau_norm = (tau_smooth - np.mean(tau_smooth)) / max(1e-4, np.std(tau_smooth))
    fossil_fluct = np.clip(tau_norm * 0.4, -0.8, 0.8)
    
    # 2. Nuevas fluctuaciones cuánticas primordiales
    noise_fft = np.fft.fft2(np.random.randn(GRID_SIZE, GRID_SIZE))
    fluct_new = np.real(np.fft.ifft2(noise_fft * p_k))
    fluct_new = (fluct_new - np.mean(fluct_new)) / np.std(fluct_new) * 0.35
    
    # 3. Plasma expulsado del Agujero Blanco + cúmulos galácticos secundarios
    core_blast = 1.8 * np.exp(-(r**2) / 35.0)
    shell_blast = 1.5 * np.exp(-((r - 12.0)**2) / 25.0)
    
    # Cúmulo galáctico secundario fragmentado por la onda de choque
    cx2, cy2 = (x0 + np.random.randint(15, 30)) % GRID_SIZE, (y0 + np.random.randint(15, 30)) % GRID_SIZE
    dx2 = (X - cx2 + GRID_SIZE / 2) % GRID_SIZE - GRID_SIZE / 2
    dy2 = (Y - cy2 + GRID_SIZE / 2) % GRID_SIZE - GRID_SIZE / 2
    r2 = np.sqrt(dx2**2 + dy2**2)
    seed2 = 1.6 * np.exp(-(r2**2) / 16.0)
    
    rho_new = np.maximum(0.06, 0.9 + fluct_new + fossil_fluct + core_blast + shell_blast + seed2)
    
    # 4. Impulso Cinético Radial de Expansión (limitado estrictamente por c)
    v_blast = np.minimum(C_LIGHT * 0.75, 1.8 * np.exp(-r / 15.0))
    vx_blast = v_blast * (dx / r_safe)
    vy_blast = v_blast * (dy / r_safe)
    
    # 5. Recalentamiento Térmico Primordial (Alta presión de radiación post-inflación)
    T_new = 25.0 * (rho_new**0.5) + 35.0
    
    return rho_new, vx_blast, vy_blast, T_new

# =====================================================================
# MOTOR FÍSICO COSMOLÓGICO CON REOTRANSDUCTOR, INFLACIÓN Y RELATIVIDAD
# =====================================================================
def update_cosmology(rho, T, I, tau, tau_eon_start, v_x, v_y, scale_factor, eon):
    # a. Expansión cósmica con Fase de Inflación Primordial (Super-expansión del espacio)
    if scale_factor < 1.05:
        # Tasa de Hubble acelerada por repulsión cuántica de Planck / Inflatón
        H_eff = H_0 * (1.0 + INFLATION_BOOST * np.exp(-(scale_factor - 1.0) / 0.015))
    else:
        H_eff = H_0
        
    scale_factor += H_eff * DT
    
    # b. Gravedad: Ecuación de Poisson ∇²Φ = 4πG(ρ - ρ̄) vía FFT 2D
    delta_rho = rho - np.mean(rho)
    delta_rho_fft = np.fft.fft2(delta_rho)
    phi_fft = -4.0 * np.pi * G_CONST * delta_rho_fft / k2
    phi_fft[0, 0] = 0.0
    phi = np.real(np.fft.ifft2(phi_fft))
    
    grad_phi_y, grad_phi_x = np.gradient(phi)
    
    # c. Gradiente de Presión de Jeans (Polítropo P = CS2 * ρ^1.3)
    P = CS2 * (rho**1.3)
    grad_P_y, grad_P_x = np.gradient(P)
    
    # Aceleración neta: Gravedad - Gradiente de presión
    acc_x = -grad_phi_x - (grad_P_x / rho)
    acc_y = -grad_phi_y - (grad_P_y / rho)
    
    # Evolución hidrodinámica con inercia cinética y amortiguamiento de Hubble
    v_x = 0.92 * v_x + 0.08 * (0.06 * acc_x)
    v_y = 0.92 * v_y + 0.08 * (0.06 * acc_y)
    
    # Límite Causal Relativista de la Velocidad de la Materia (v <= c)
    v_mag = np.sqrt(v_x**2 + v_y**2)
    v_limit = np.maximum(1.0, v_mag / C_LIGHT)
    v_x = v_x / v_limit
    v_y = v_y / v_limit
    
    # d. Hidrodinámica de materia cósmica
    flux_x = rho * v_x
    flux_y = rho * v_y
    div_flux = np.gradient(flux_x, axis=1) + np.gradient(flux_y, axis=0)
    
    laplacian_rho = (
        np.roll(rho, 1, axis=0) + np.roll(rho, -1, axis=0) +
        np.roll(rho, 1, axis=1) + np.roll(rho, -1, axis=1) - 4 * rho
    )
    rho_next = np.clip(rho - div_flux * DT + 0.04 * laplacian_rho * DT, 0.02, 12.0)
    
    # e. Evolución térmica: Difusión + Compresión - Enfriamiento de Hubble
    laplacian_T = (
        np.roll(T, 1, axis=0) + np.roll(T, -1, axis=0) +
        np.roll(T, 1, axis=1) + np.roll(T, -1, axis=1) - 4 * T
    )
    compression_heating = 6.0 * np.maximum(0.0, -div_flux)
    hubble_cooling = H_eff * T
    T_next = np.clip(T + (DIFFUSION_COEFF * laplacian_T + compression_heating - hubble_cooling) * DT, 2.73, 2000.0)
    
    # f. Corrientes Disipativas y Producción de Entropía de Onsager
    inv_T = 1.0 / T_next
    grad_inv_T_y, grad_inv_T_x = np.gradient(inv_T)
    grad_T_y, grad_T_x = np.gradient(T_next)
    
    J_T_x = -DIFFUSION_COEFF * grad_T_x
    J_T_y = -DIFFUSION_COEFF * grad_T_y
    
    sigma_thermal = np.maximum(0.0, J_T_x * grad_inv_T_x + J_T_y * grad_inv_T_y)
    sigma_grav = (rho_next * (grad_phi_x**2 + grad_phi_y**2)) / (T_next * 50.0)
    sigma_total = sigma_thermal + sigma_grav
    
    # g. Reotransductor: El tiempo propio acumulado nunca se borra (acumulador continuo eterno)
    d_tau_dt = KAPPA * sigma_total
    tau_next = tau + d_tau_dt * DT
    
    # Tiempo y entropía acumulados específicamente en el EÓN ACTUAL
    tau_current_eon = tau_next - tau_eon_start
    
    # h. Campo Informacional Continuo I(r, t)
    flux_I_x = I * v_x
    flux_I_y = I * v_y
    div_flux_I = np.gradient(flux_I_x, axis=1) + np.gradient(flux_I_y, axis=0)
    
    laplacian_I = (
        np.roll(I, 1, axis=0) + np.roll(I, -1, axis=0) +
        np.roll(I, 1, axis=1) + np.roll(I, -1, axis=1) - 4 * I
    )
    
    sustenance = 0.6 * sigma_total * (rho_next / np.mean(rho_next))
    thermal_noise = 0.0004 * T_next
    dI_dt = -div_flux_I + 0.02 * laplacian_I + (sustenance - thermal_noise - LANDAUER_DECAY * I)
    I_next = np.clip(I + dI_dt * DT, 0.0, 1.0)
    
    # i. CONDICIÓN FÍSICA DE REBOTE POR LÍMITE DE BEKENSTEIN EN EL EÓN ACTUAL
    total_mass = np.sum(rho_next)
    core_mask = rho_next > 1.0
    core_mass = np.sum(rho_next[core_mask])
    mass_fraction = core_mass / total_mass
    
    s_bh_eon = np.max(tau_current_eon)
    tau_bekenstein_crit = ZETA_BEKENSTEIN * max(1.0, (core_mass / 30.0)**2)
    tunnel_progress = min(1.0, s_bh_eon / max(1.0, tau_bekenstein_crit))
    
    # Disparo del Agujero Blanco cuando el colapso gravitatorio es maduro (35%+ de masa concentrada)
    # y la entropía acumulada DURANTE ESTE EÓN satura el límite de Bekenstein
    if mass_fraction >= 0.35 and s_bh_eon >= tau_bekenstein_crit:
        eon += 1
        scale_factor = 1.0
        tau_eon_start = tau_next.copy()  # Establece la nueva línea base para el nuevo eón
        rho_next, v_x, v_y, T_next = trigger_white_hole_eon(rho_next, tau_next)
        I_next = np.clip((rho_next - 0.5) / 2.5, 0.0, 1.0)
    
    r_hubble = C_LIGHT / max(1e-5, H_eff)
    return rho_next, T_next, I_next, tau_next, tau_eon_start, v_x, v_y, scale_factor, eon, d_tau_dt, phi, mass_fraction, s_bh_eon, tau_bekenstein_crit, tunnel_progress, r_hubble

# =====================================================================
# INTERFAZ GRÁFICA Y DASHBOARD DE 5 PANELES + CONSOLA DE TELEMETRÍA
# =====================================================================
fig = plt.figure(figsize=(14, 8.8))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.28, wspace=0.30, top=0.92, bottom=0.08, left=0.05, right=0.95)

ax_rho     = fig.add_subplot(gs[0, 0])
ax_rate    = fig.add_subplot(gs[0, 1])
ax_index   = fig.add_subplot(gs[0, 2])
ax_tau     = fig.add_subplot(gs[1, 0])
ax_log_tau = fig.add_subplot(gs[1, 1])
ax_info    = fig.add_subplot(gs[1, 2])

# Panel 1: Densidad Cósmica
im_rho = ax_rho.imshow(rho, cmap='magma', origin='lower', vmin=0.0, vmax=5.0)
title_rho = ax_rho.set_title("1. Densidad de Materia Cósmica (ρ)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_rho, ax=ax_rho, label="Densidad Normalizada [ρ/ρ̄]", fraction=0.046, pad=0.04)

# Panel 2: Velocidad del Reotransductor
im_rate = ax_rate.imshow(np.zeros((GRID_SIZE, GRID_SIZE)), cmap='plasma', origin='lower', vmin=0.0, vmax=2.0)
ax_rate.set_title("2. Velocidad Reotransductor (dτ/dt)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_rate, ax=ax_rate, label="Flujo Local de Presente", fraction=0.046, pad=0.04)

# Panel 3: Autoorganización / Negentropía Continua
im_index = ax_index.imshow(I, cmap='cividis', origin='lower', vmin=0.0, vmax=1.0)
title_index = ax_index.set_title("3. Autoorganización I(r,t)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_index, ax=ax_index, label="Orden Informacional [0 a 1]", fraction=0.046, pad=0.04)

# Panel 4: Tiempo Propio Emergente Lineal
im_tau = ax_tau.imshow(tau, cmap='viridis', origin='lower', vmin=0.0, vmax=50.0)
title_tau = ax_tau.set_title("4. Tiempo Emergente Lineal (τ)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_tau, ax=ax_tau, label="Segundos Cosmológicos τ", fraction=0.046, pad=0.04)

# Panel 5: Red Fósil Logarítmica (Arqueología Cósmica)
im_log_tau = ax_log_tau.imshow(np.log10(1.0 + tau), cmap='inferno', origin='lower', vmin=0.0, vmax=5.0)
title_log_tau = ax_log_tau.set_title("5. Memoria Fósil Logarítmica log₁₀(1+τ)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_log_tau, ax=ax_log_tau, label="Escala Logarítmica de Tiempo", fraction=0.046, pad=0.04)

# Panel 6: Consola de Telemetría de Bekenstein & Eones
ax_info.axis('off')
title_info = ax_info.text(0.5, 0.96, "TELEMETRÍA DE BEKENSTEIN & EONES", fontweight='bold', fontsize=10.0, ha='center', va='top', color='#0f172a')
telemetry_box = ax_info.text(
    0.04, 0.82, "", fontsize=8.2, fontfamily='monospace', va='top',
    bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8fafc', edgecolor='#cbd5e1')
)

# Slider de Aceleración (hasta 1000x)
ax_slider = fig.add_axes([0.25, 0.02, 0.50, 0.025])
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

# Bucle de animación continuo
def animate(frame):
    global rho, T, I, tau, tau_eon_start, v_x, v_y, scale_factor, eon
    
    mass_frac = 0.0
    s_bh_val = 0.0
    s_crit = 0.0
    progress = 0.0
    r_hub = 0.0
    
    for _ in range(steps_per_frame):
        rho, T, I, tau, tau_eon_start, v_x, v_y, scale_factor, eon, d_tau_dt, phi, mass_frac, s_bh_val, s_crit, progress, r_hub = update_cosmology(
            rho, T, I, tau, tau_eon_start, v_x, v_y, scale_factor, eon
        )
    
    # Actualizar matrices gráficas
    im_rho.set_array(rho)
    im_rate.set_array(d_tau_dt)
    im_index.set_array(I)
    im_tau.set_array(tau)
    im_log_tau.set_array(np.log10(1.0 + tau))
    
    # Identificación de la era interna del Eón
    if scale_factor < 1.05:
        era_str = "Fase de Inflación Cuántica Primordial"
    elif scale_factor < 2.5:
        era_str = "Era de Galaxias y Filamentos"
    elif scale_factor < 7.0:
        era_str = "Era de Fusiones y Maduración"
    elif mass_frac >= 0.35:
        era_str = "Era del Agujero Negro Virializado"
    else:
        era_str = "Fase Asintótica Pre-Rebote"
        
    redshift = max(0.0, (4.50 / max(0.01, float(scale_factor))) - 1.0)
    fig.suptitle(f"Simulación Cosmológica: Reotransductor, Bekenstein & Agujero Blanco | Eón N = {eon} [{era_str}]", fontsize=11.5, fontweight='bold')
    
    title_rho.set_text(f"1. Densidad Materia (ρ) [a = {scale_factor:.2f}]")
    title_index.set_text(f"3. Autoorganización I(r,t) [I_max: {np.max(I):.2f}]")
    title_tau.set_text(f"4. Tiempo Lineal (τ) [Pico: {np.max(tau):.0f} s]")
    title_log_tau.set_text(f"5. Red Fósil log₁₀(1+τ) [Max: {np.max(np.log10(1.0 + tau)):.2f}]")
    
    # Telemetría en consola interactiva
    bar_len = 10
    filled_len = int(progress * bar_len)
    bar_str = "█" * filled_len + "░" * (bar_len - filled_len)
    
    tau_core_val = np.max(tau)
    tau_void_val = np.mean(tau[36:40, 8:12])
    desynch_ratio = tau_core_val / max(1e-4, tau_void_val)
    
    telemetry_str = (
        f"• Eón Cósmico Actual: N = {eon}\n"
        f"• Factor de Escala: a = {scale_factor:.3f} (z = {redshift:.2f})\n"
        f"• Velocidad de la Luz (c): {C_LIGHT:.2f} celdas/s\n"
        f"• Masa Núcleo (M_BH): {mass_frac * 100.0:.1f}% del total\n"
        f"• Entropía Eón N (S_BH): {s_bh_val:.0f} k_B\n"
        f"• Límite Bekenstein (S_max): {s_crit:.0f} k_B\n"
        f"• Túnel Cuántico Eón {eon}: [{bar_str}] {progress * 100.0:.1f}%\n"
        f"• Desincronización Cúmulo/Vacío: {desynch_ratio:.1f}x\n"
        f"• Odómetro Fósil Acumulado: {np.max(tau):.0f} s\n"
        f"• Estado: {'Agujero Blanco Inminente' if progress >= 0.95 else 'Evolución Hidrodinámica'}"
    )
    telemetry_box.set_text(telemetry_str)
    
    # Escalas dinámicas
    im_rho.set_clim(vmin=0.0, vmax=max(3.0, np.percentile(rho, 99)))
    im_rate.set_clim(vmin=0.0, vmax=max(0.5, np.percentile(d_tau_dt, 99)))
    im_tau.set_clim(vmin=0.0, vmax=max(1.0, np.max(tau)))
    im_log_tau.set_clim(vmin=0.0, vmax=max(1.0, np.max(np.log10(1.0 + tau))))
    
    return im_rho, im_rate, im_index, im_tau, im_log_tau, telemetry_box

if __name__ == '__main__':
    ani = FuncAnimation(fig, animate, interval=30, blit=False, cache_frame_data=False)
    plt.show()
