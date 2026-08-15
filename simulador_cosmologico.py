import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
import matplotlib.patches as mpatches

# =====================================================================
# CONFIGURACIÓN DEL MODELO COSMOLÓGICO Y PARÁMETROS FÍSICOS
# =====================================================================
GRID_SIZE = 50              # Resolución de la malla espacial cosmológica (N x N)
DT = 0.05                   # Paso temporal de integración
DIFFUSION_COEFF = 0.3       # Conductividad térmica del plasma intergaláctico (k)
KAPPA = 50.0                # Constante de acoplamiento del Reotransductor (κ)
LANDAUER_DECAY = 0.015      # Tasa de decaimiento entrópico de estructuras complejas (γ)
G_CONST = 0.05              # Constante de gravitación efectiva normalizada (G)
H_0 = 0.0003                # Tasa de expansión de Hubble inicial

# 1. Campo de densidad primordial con espectro de potencias de perturbaciones P(k) ~ k^(-0.7)
np.random.seed(137)
kx = 2 * np.pi * np.fft.fftfreq(GRID_SIZE)[:, None]
ky = 2 * np.pi * np.fft.fftfreq(GRID_SIZE)[None, :]
k2 = kx**2 + ky**2
k2[0, 0] = 1.0  # Evitar división por cero en modo k=0

noise_fft = np.fft.fft2(np.random.randn(GRID_SIZE, GRID_SIZE))
p_k = 1.0 / (k2**0.7)
p_k[0, 0] = 0.0
fluctuations = np.real(np.fft.ifft2(noise_fft * p_k))
fluctuations = (fluctuations - np.mean(fluctuations)) / np.std(fluctuations) * 0.5

# Densidad de materia cósmica ρ(x, y)
rho = np.maximum(0.1, 1.0 + fluctuations)
rho[22:27, 22:27] += 3.5    # Semilla A: Supercúmulo central denso
rho[36:40, 12:16] += 2.2    # Semilla B: Filamento galáctico activo
rho[8:12, 38:42] = 0.08     # Región C: Vacío cósmico profundo (underdense void)

# 2. Campo de Temperatura T(x, y) acoplado a la densidad (bariones y radiación)
T = 15.0 * (rho**0.5) + 2.73

# 3. Índices de Negentropía Cósmica I(x, y) (Estructuras de orden, ej: galaxias y biosferas)
I = np.zeros((GRID_SIZE, GRID_SIZE))
I[22:26, 22:26] = 1.0       # Estructura A: Supercúmulo
I[36:40, 12:16] = 1.0       # Estructura B: Galaxia en filamento
I[8:12, 38:42] = 1.0        # Estructura C: Galaxia enana aislada en el vacío

# 4. Tiempo Térmico Propio Emergente τ(x, y) y Factor de Escala a(t)
tau = np.zeros((GRID_SIZE, GRID_SIZE))
scale_factor = 1.0
steps_per_frame = 2

# =====================================================================
# MOTOR FÍSICO COSMOLÓGICO CON REOTRANSDUCTOR LOCAL
# =====================================================================
def update_cosmology(rho, T, I, tau, scale_factor):
    # a. Expansión cosmológica
    scale_factor += H_0 * DT
    
    # b. Gravedad: Solución exacta de la ecuación de Poisson ∇²Φ = 4πG(ρ - ρ̄) vía FFT 2D
    delta_rho = rho - np.mean(rho)
    delta_rho_fft = np.fft.fft2(delta_rho)
    phi_fft = -4.0 * np.pi * G_CONST * delta_rho_fft / k2
    phi_fft[0, 0] = 0.0
    phi = np.real(np.fft.ifft2(phi_fft))
    
    # Aceleración gravitatoria: g = -∇Φ
    grad_phi_y, grad_phi_x = np.gradient(phi)
    v_x = -0.05 * np.clip(grad_phi_x, -5.0, 5.0)
    v_y = -0.05 * np.clip(grad_phi_y, -5.0, 5.0)
    
    # c. Hidrodinámica de materia cósmica: Ecuación de continuidad + difusión Jeans
    flux_x = rho * v_x
    flux_y = rho * v_y
    div_flux = np.gradient(flux_x, axis=1) + np.gradient(flux_y, axis=0)
    
    laplacian_rho = (
        np.roll(rho, 1, axis=0) + np.roll(rho, -1, axis=0) +
        np.roll(rho, 1, axis=1) + np.roll(rho, -1, axis=1) - 4 * rho
    )
    rho_next = np.clip(rho - div_flux * DT + 0.05 * laplacian_rho * DT, 0.02, 25.0)
    
    # d. Evolución térmica: Conducción Fourier + Calentamiento adiabático por compresión - Enfriamiento de Hubble
    laplacian_T = (
        np.roll(T, 1, axis=0) + np.roll(T, -1, axis=0) +
        np.roll(T, 1, axis=1) + np.roll(T, -1, axis=1) - 4 * T
    )
    compression_heating = 8.0 * np.maximum(0.0, -div_flux)
    hubble_cooling = H_0 * T
    T_next = np.clip(T + (DIFFUSION_COEFF * laplacian_T + compression_heating - hubble_cooling) * DT, 2.73, 3000.0)
    
    # e. Corrientes Disipativas y Producción de Entropía de Onsager
    inv_T = 1.0 / T_next
    grad_inv_T_y, grad_inv_T_x = np.gradient(inv_T)
    grad_T_y, grad_T_x = np.gradient(T_next)
    
    J_T_x = -DIFFUSION_COEFF * grad_T_x
    J_T_y = -DIFFUSION_COEFF * grad_T_y
    
    # Producción local de entropía térmica + gravitatoria
    sigma_thermal = np.maximum(0.0, J_T_x * grad_inv_T_x + J_T_y * grad_inv_T_y)
    sigma_grav = (rho_next * (grad_phi_x**2 + grad_phi_y**2)) / (T_next * 80.0)
    sigma_total = sigma_thermal + sigma_grav
    
    # f. Reotransductor Cosmológico: Velocidad de emergencia del tiempo propio local
    d_tau_dt = KAPPA * sigma_total
    tau_next = tau + d_tau_dt * DT
    
    # g. Dinámica de Negentropía Estructural (Límite de Landauer)
    # Las galaxias se sostienen con el flujo libre disipado de acreción y calor
    energy_flux = np.sqrt(J_T_x**2 + J_T_y**2) + 0.3 * np.sqrt(flux_x**2 + flux_y**2)
    thermal_noise = 0.001 * T_next
    sustenance = 4.5 * energy_flux
    
    mask = I > 0
    delta_I = (sustenance - thermal_noise - LANDAUER_DECAY) * DT
    I_next = np.clip(I + delta_I * mask, 0.0, 1.0)
    
    return rho_next, T_next, I_next, tau_next, scale_factor, d_tau_dt, phi

# =====================================================================
# INTERFAZ GRÁFICA Y DASHBOARD COSMOLÓGICO
# =====================================================================
fig, axes = plt.subplots(2, 2, figsize=(10, 9.2))
fig.suptitle("Simulación Cosmológica: Reotransductor del Presente Activo", fontsize=13, fontweight='bold')

ax_rho   = axes[0, 0]
ax_rate  = axes[0, 1]
ax_tau   = axes[1, 0]
ax_index = axes[1, 1]

# Panel 1: Densidad Cósmica
im_rho = ax_rho.imshow(rho, cmap='magma', origin='lower', vmin=0.0, vmax=5.0)
title_rho = ax_rho.set_title("Densidad de Materia Cósmica (ρ) [a=1.00]", fontweight='bold', fontsize=9.5)
fig.colorbar(im_rho, ax=ax_rho, label="Densidad Normalizada [ρ/ρ̄]", fraction=0.046, pad=0.04)

# Panel 2: Velocidad del Reotransductor
im_rate = ax_rate.imshow(np.zeros((GRID_SIZE, GRID_SIZE)), cmap='plasma', origin='lower', vmin=0.0, vmax=2.0)
ax_rate.set_title("Velocidad del Reotransductor (dτ/dt)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_rate, ax=ax_rate, label="Flujo Local de Presente", fraction=0.046, pad=0.04)

# Panel 3: Tiempo Propio Emergente Acumulado
im_tau = ax_tau.imshow(tau, cmap='viridis', origin='lower', vmin=0.0, vmax=50.0)
title_tau = ax_tau.set_title("Tiempo Propio Emergente (Coordenada τ)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_tau, ax=ax_tau, label="Segundos Cosmológicos Emergentes", fraction=0.046, pad=0.04)

# Panel 4: Integridad de Negentropía Cósmica
im_index = ax_index.imshow(I, cmap='cividis', origin='lower', vmin=0.0, vmax=1.0)
ax_index.set_title("Integridad Negentrópica Cósmica (I)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_index, ax=ax_index, label="Orden Estructural [0 a 1]", fraction=0.046, pad=0.04)

# Recuadros y etiquetas para las estructuras cósmicas A, B y C
ax_index.add_patch(mpatches.Rectangle((21.5, 21.5), 4, 4, linewidth=1.5, edgecolor='#00ffff', facecolor='none'))
ax_index.add_patch(mpatches.Rectangle((11.5, 35.5), 4, 4, linewidth=1.5, edgecolor='#39ff14', facecolor='none'))
ax_index.add_patch(mpatches.Rectangle((37.5, 7.5), 4, 4, linewidth=1.5, edgecolor='#ff3131', facecolor='none'))

ax_index.text(23.5, 26.5, 'A', color='#00ffff', fontweight='bold', fontsize=10, ha='center', va='bottom')
ax_index.text(13.5, 40.5, 'B', color='#39ff14', fontweight='bold', fontsize=10, ha='center', va='bottom')
ax_index.text(39.5, 12.5, 'C', color='#ff3131', fontweight='bold', fontsize=10, ha='center', va='bottom')

# Leyenda explicativa
legend_elements = [
    mpatches.Patch(facecolor='none', edgecolor='#00ffff', linewidth=1.5, label='A: Supercúmulo Central - Alto colapso y flujo negentrópico'),
    mpatches.Patch(facecolor='none', edgecolor='#39ff14', linewidth=1.5, label='B: Filamento Galáctico - Sostenido por acreción moderada'),
    mpatches.Patch(facecolor='none', edgecolor='#ff3131', linewidth=1.5, label='C: Galaxia en Vacío Cósmico - Decae por Landauer')
]
ax_index.legend(handles=legend_elements, loc='upper left', fontsize=7.2, framealpha=0.85)

# Ajuste fino de la cuadrícula
plt.subplots_adjust(top=0.93, bottom=0.09, left=0.06, right=0.94, hspace=0.22, wspace=0.28)

# Slider de Aceleración en la parte inferior
ax_slider = fig.add_axes([0.25, 0.025, 0.62, 0.03])
slider_speed = Slider(
    ax=ax_slider,
    label='Aceleración (Pasos / Frame) ',
    valmin=1,
    valmax=50,
    valinit=2,
    valstep=1,
    color='#2563eb'
)

def update_speed(val):
    global steps_per_frame
    steps_per_frame = int(slider_speed.val)

slider_speed.on_changed(update_speed)

# Bucle de animación
def animate(frame):
    global rho, T, I, tau, scale_factor
    
    for _ in range(steps_per_frame):
        rho, T, I, tau, scale_factor, d_tau_dt, phi = update_cosmology(rho, T, I, tau, scale_factor)
    
    # Actualizar matrices gráficas
    im_rho.set_array(rho)
    im_rate.set_array(d_tau_dt)
    im_tau.set_array(tau)
    im_index.set_array(I)
    
    # Redshift z = 1/a - 1
    redshift = max(0.0, (1.0 / scale_factor) - 1.0)
    title_rho.set_text(f"Densidad Materia (ρ) | a = {scale_factor:.3f} (z = {redshift:.2f})")
    
    # Telemetría de desincronización de tiempo cósmico
    tau_void = np.mean(tau[8:12, 38:42])
    tau_cluster = np.mean(tau[22:26, 22:26])
    ratio = (tau_cluster / max(1e-4, tau_void))
    title_tau.set_text(f"Tiempo Emergente (τ) [Cúmulo/Vacío: {ratio:.1f}x]")
    
    # Escalas dinámicas
    im_rho.set_clim(vmin=0.0, vmax=max(3.0, np.percentile(rho, 99)))
    im_rate.set_clim(vmin=0.0, vmax=max(0.5, np.percentile(d_tau_dt, 99)))
    im_tau.set_clim(vmin=0.0, vmax=max(1.0, np.max(tau)))
    
    return im_rho, im_rate, im_tau, im_index, title_rho, title_tau

ani = FuncAnimation(fig, animate, interval=30, blit=False, cache_frame_data=False)
plt.show()
