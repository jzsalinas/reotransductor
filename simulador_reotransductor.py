import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
import matplotlib.patches as mpatches

# =====================================================================
# CONFIGURACIÓN DE PARÁMETROS FÍSICOS
# =====================================================================
GRID_SIZE = 50          # Resolución espacial
DT = 0.1                # Intervalo temporal estable (CFL)
DIFFUSION_COEFF = 0.5   # Conductividad térmica (k)
KAPPA = 50.0            # Constante disipativa de Onsager
LANDAUER_DECAY = 0.02   # Tasa de decaimiento entrópico de los Índices

# Reservas iniciales de combustible térmico (Joules de la caldera)
FUEL_1_INIT = 35000.0   # Caldera 1 (1000 K)
FUEL_2_INIT = 25000.0   # Caldera 2 (800 K)

fuel_1 = FUEL_1_INIT
fuel_2 = FUEL_2_INIT

# 1. Temperatura inicial (Vacío frío a 2.73 K)
T = np.ones((GRID_SIZE, GRID_SIZE)) * 2.73
T[15:18, 15:18] = 1000.0
T[32:35, 32:35] = 800.0

# 2. Índices de Negentropía iniciales
I = np.zeros((GRID_SIZE, GRID_SIZE))
I[18:22, 18:22] = 1.0  # Estructura A
I[35:39, 35:39] = 1.0  # Estructura B
I[5:9, 40:44] = 1.0    # Estructura C (Periferia)

# 3. Tiempo Térmico Acumulado
tau = np.zeros((GRID_SIZE, GRID_SIZE))
steps_per_frame = 1

# =====================================================================
# MOTOR FÍSICO CON RECURSOS FINITOS
# =====================================================================
def update_physics(T, I, tau):
    global fuel_1, fuel_2
    
    # a. Conducción de Fourier (Diferencias finitas en volumen interno)
    laplacian = np.zeros_like(T)
    laplacian[1:-1, 1:-1] = (
        T[:-2, 1:-1] + T[2:, 1:-1] +
        T[1:-1, :-2] + T[1:-1, 2:] - 4 * T[1:-1, 1:-1]
    )
    T_next = T + DIFFUSION_COEFF * laplacian * DT
    
    # b. Inyección y enfriamiento gradual según combustible (Reservorio Lineal)
    # T_caldera(t) = T_base + (T_max - T_base) * (combustible / combustible_inicial)
    pct_f1 = fuel_1 / FUEL_1_INIT
    target_T1 = 2.73 + (1000.0 - 2.73) * pct_f1
    if fuel_1 > 0:
        burn_rate_1 = np.sum(np.maximum(0.0, target_T1 - T_next[15:18, 15:18])) * 0.5 * DT
        fuel_1 = max(0.0, fuel_1 - burn_rate_1)
        pct_f1 = fuel_1 / FUEL_1_INIT
        T_next[15:18, 15:18] = 2.73 + (1000.0 - 2.73) * pct_f1
    
    pct_f2 = fuel_2 / FUEL_2_INIT
    target_T2 = 2.73 + (800.0 - 2.73) * pct_f2
    if fuel_2 > 0:
        burn_rate_2 = np.sum(np.maximum(0.0, target_T2 - T_next[32:35, 32:35])) * 0.5 * DT
        fuel_2 = max(0.0, fuel_2 - burn_rate_2)
        pct_f2 = fuel_2 / FUEL_2_INIT
        T_next[32:35, 32:35] = 2.73 + (800.0 - 2.73) * pct_f2
        
    # Fronteras frías fijas (Sumidero cósmico a 2.73 K)
    T_next[0, :] = T_next[-1, :] = T_next[:, 0] = T_next[:, -1] = 2.73
    
    # c. Gradiente inverso y corriente disipativa
    inv_T = 1.0 / T_next
    grad_inv_T_y, grad_inv_T_x = np.gradient(inv_T)
    
    grad_T_y, grad_T_x = np.gradient(T_next)
    J_x = -DIFFUSION_COEFF * grad_T_x
    J_y = -DIFFUSION_COEFF * grad_T_y
    
    # d. Velocidad del Reotransductor (Producción de Entropía)
    d_tau_dt = KAPPA * np.maximum(0.0, (J_x * grad_inv_T_x + J_y * grad_inv_T_y))
    tau_next = tau + d_tau_dt * DT
    
    # e. Dinámica de Supervivencia de los Índices
    energy_flux = np.sqrt(J_x**2 + J_y**2)
    thermal_noise = 0.001 * T_next
    sustenance = 5.0 * energy_flux
    
    # Solo las estructuras no extintas (I > 0) evolucionan
    mask = I > 0
    delta_I = np.zeros_like(I)
    delta_I[mask] = (sustenance[mask] - thermal_noise[mask] - LANDAUER_DECAY) * DT
    I_next = np.clip(I + delta_I, 0.0, 1.0)
    
    return T_next, I_next, tau_next, d_tau_dt

# =====================================================================
# INTERFAZ GRÁFICA Y DASHBOARD
# =====================================================================
fig, axes = plt.subplots(2, 2, figsize=(9.5, 9.2))
fig.suptitle("Simulación Termodinámica: Recursos Finitos y Muerte del Tiempo", fontsize=13, fontweight='bold')

ax_temp  = axes[0, 0]
ax_rate  = axes[0, 1]
ax_tau   = axes[1, 0]
ax_index = axes[1, 1]

# Plot 1: Temperatura
im_temp = ax_temp.imshow(T, cmap='inferno', origin='lower', vmin=2.73, vmax=1000.0)
title_temp = ax_temp.set_title(f"Temperatura (C1: 100% | C2: 100%)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_temp, ax=ax_temp, label="Kelvins [K]", fraction=0.046, pad=0.04)

# Plot 2: Velocidad del Reotransductor (d_tau/dt)
im_rate = ax_rate.imshow(np.zeros((GRID_SIZE, GRID_SIZE)), cmap='plasma', origin='lower', vmin=0, vmax=8.0)
ax_rate.set_title(r"Velocidad del Reotransductor ($d\tau/dt$)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_rate, ax=ax_rate, label="Flujo de Presente", fraction=0.046, pad=0.04)

# Plot 3: Tiempo Emergente Acumulado (tau)
im_tau = ax_tau.imshow(tau, cmap='viridis', origin='lower', vmin=0, vmax=100.0)
ax_tau.set_title(r"Tiempo Acumulado (Coordenada $\tau$)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_tau, ax=ax_tau, label="Segundos Emergentes", fraction=0.046, pad=0.04)

# Plot 4: Integridad de los Índices de Negentropía (I)
im_index = ax_index.imshow(I, cmap='cividis', origin='lower', vmin=0.0, vmax=1.0)
ax_index.set_title(r"Integridad del Índice ($I$)", fontweight='bold', fontsize=9.5)
fig.colorbar(im_index, ax=ax_index, label="Orden Estructural [0 a 1]", fraction=0.046, pad=0.04)

# Recuadros y etiquetas para las estructuras A, B y C
ax_index.add_patch(mpatches.Rectangle((17.5, 17.5), 4, 4, linewidth=1.5, edgecolor='#00ffff', facecolor='none'))
ax_index.add_patch(mpatches.Rectangle((34.5, 34.5), 4, 4, linewidth=1.5, edgecolor='#39ff14', facecolor='none'))
ax_index.add_patch(mpatches.Rectangle((39.5, 4.5), 4, 4, linewidth=1.5, edgecolor='#ff3131', facecolor='none'))

ax_index.text(19.5, 22.5, 'A', color='#00ffff', fontweight='bold', fontsize=10, ha='center', va='bottom')
ax_index.text(36.5, 39.5, 'B', color='#39ff14', fontweight='bold', fontsize=10, ha='center', va='bottom')
ax_index.text(41.5, 9.5, 'C', color='#ff3131', fontweight='bold', fontsize=10, ha='center', va='bottom')

# Leyenda explicativa de los índices
legend_elements = [
    mpatches.Patch(facecolor='none', edgecolor='#00ffff', linewidth=1.5, label='A: C1 (1000 K) - Sostenida por alto flujo'),
    mpatches.Patch(facecolor='none', edgecolor='#39ff14', linewidth=1.5, label='B: C2 (800 K) - Sostenida por flujo medio'),
    mpatches.Patch(facecolor='none', edgecolor='#ff3131', linewidth=1.5, label='C: Vacío (2.73 K) - Decae por Landauer')
]
ax_index.legend(handles=legend_elements, loc='upper left', fontsize=7.2, framealpha=0.85)

# Ajuste fino de la cuadrícula de subplots
plt.subplots_adjust(top=0.93, bottom=0.09, left=0.06, right=0.94, hspace=0.22, wspace=0.28)

# Slider de Aceleración en la parte inferior
ax_slider = fig.add_axes([0.25, 0.025, 0.62, 0.03])
slider_speed = Slider(
    ax=ax_slider,
    label='Aceleración (Pasos / Frame) ',
    valmin=1,
    valmax=50,
    valinit=5,
    valstep=1,
    color='#2563eb'
)

def update_speed(val):
    global steps_per_frame
    steps_per_frame = int(slider_speed.val)

slider_speed.on_changed(update_speed)

# Bucle de animación
def animate(frame):
    global T, I, tau
    
    for _ in range(steps_per_frame):
        T, I, tau, d_tau_dt = update_physics(T, I, tau)
    
    # Actualizar matrices gráficas
    im_temp.set_array(T)
    im_rate.set_array(d_tau_dt)
    im_tau.set_array(tau)
    im_index.set_array(I)
    
    # Actualizar porcentaje de combustible y temperatura actual en el título
    pct_f1 = (fuel_1 / FUEL_1_INIT) * 100.0
    pct_f2 = (fuel_2 / FUEL_2_INIT) * 100.0
    t_c1 = 2.73 + (1000.0 - 2.73) * (fuel_1 / FUEL_1_INIT)
    t_c2 = 2.73 + (800.0 - 2.73) * (fuel_2 / FUEL_2_INIT)
    title_temp.set_text(f"Temperatura (C1: {pct_f1:.1f}% [{t_c1:.0f} K] | C2: {pct_f2:.1f}% [{t_c2:.0f} K])")
    
    # Escala de tiempo dinámico
    im_tau.set_clim(vmin=0, vmax=max(1.0, np.max(tau)))
    
    return im_temp, im_rate, im_tau, im_index, title_temp

ani = FuncAnimation(fig, animate, interval=30, blit=False, cache_frame_data=False)
plt.show()