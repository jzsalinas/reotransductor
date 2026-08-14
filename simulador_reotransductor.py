import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as mpatches

# =====================================================================
# CONFIGURACIÓN DE PARÁMETROS FÍSICOS Y GEOMETRÍA
# =====================================================================
GRID_SIZE = 50          # Resolución de la Malla espacial
DT = 0.1                # Intervalo de integración temporal (dt)
DIFFUSION_COEFF = 0.5   # Conductividad térmica del medio (k)
KAPPA = 50.0            # Constante disipativa del Reotransductor (kappa)
LANDAUER_DECAY = 0.02   # Tasa de decaimiento entrópico natural de los Índices

# 1. Temperatura (T): Fondo cósmico frío a 2.73 K
T = np.ones((GRID_SIZE, GRID_SIZE)) * 2.73
T[15:18, 15:18] = 1000.0  # Caldera 1
T[32:35, 32:35] = 800.0   # Caldera 2

# 2. Índices de Negentropía (I): Estructuras de orden
I = np.zeros((GRID_SIZE, GRID_SIZE))
I[18:22, 18:22] = 1.0     # Estructura A (Cerca de Caldera 1)
I[35:39, 35:39] = 1.0     # Estructura B (Cerca de Caldera 2)
I[5:9, 40:44] = 1.0       # Estructura C (Aislada en el vacío)

# 3. Coordenada de Tiempo Emergente Acumulado (tau)
tau = np.zeros((GRID_SIZE, GRID_SIZE))

# =====================================================================
# MOTOR DE EVOLUCIÓN FÍSICA (VECTORIZADO)
# =====================================================================
def update_physics(T, I, tau):
    # a. Difusión térmica con diferencias finitas sin wrap-around
    laplacian = np.zeros_like(T)
    laplacian[1:-1, 1:-1] = (
        T[:-2, 1:-1] + T[2:, 1:-1] +
        T[1:-1, :-2] + T[1:-1, 2:] - 4 * T[1:-1, 1:-1]
    )
    
    T_next = T + DIFFUSION_COEFF * laplacian * DT
    T_next[15:18, 15:18] = 1000.0
    T_next[32:35, 32:35] = 800.0
    T_next[0, :] = T_next[-1, :] = T_next[:, 0] = T_next[:, -1] = 2.73
    
    # b. Gradiente de Temperatura Inversa: X = grad(1/T)
    inv_T = 1.0 / T_next
    grad_inv_T_y, grad_inv_T_x = np.gradient(inv_T)
    
    # c. Flujo térmico: J = -k * grad(T)
    grad_T_y, grad_T_x = np.gradient(T_next)
    J_x = -DIFFUSION_COEFF * grad_T_x
    J_y = -DIFFUSION_COEFF * grad_T_y
    
    # d. Velocidad del Reotransductor (Producción de Entropía local de Onsager)
    # sigma = J . grad(1/T) >= 0
    d_tau_dt = KAPPA * (J_x * grad_inv_T_x + J_y * grad_inv_T_y)
    d_tau_dt = np.maximum(0.0, d_tau_dt)
    tau_next = tau + d_tau_dt * DT
    
    # e. Dinámica de los Índices de Negentropía
    energy_flux = np.sqrt(J_x**2 + J_y**2)
    thermal_noise = 0.001 * T_next
    sustenance = 5.0 * energy_flux
    
    mask = I > 0
    delta_I = np.zeros_like(I)
    delta_I[mask] = (sustenance[mask] - thermal_noise[mask] - LANDAUER_DECAY) * DT
    I_next = np.clip(I + delta_I, 0.0, 1.0)
    
    return T_next, I_next, tau_next, d_tau_dt

# =====================================================================
# INTERFAZ GRÁFICA
# =====================================================================
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
fig.suptitle("Simulación del Reotransductor del Presente Activo", fontsize=14, fontweight='bold')

# Plot 1: Temperatura
im_temp = axes[0, 0].imshow(T, cmap='inferno', origin='lower', vmin=2.73, vmax=1000.0)
axes[0, 0].set_title("Temperatura de la Malla ($T$)", fontweight='bold')
fig.colorbar(im_temp, ax=axes[0, 0], label="Kelvins [K]")

# Plot 2: Velocidad del Reotransductor (d_tau/dt)
im_rate = axes[0, 1].imshow(np.zeros((GRID_SIZE, GRID_SIZE)), cmap='plasma', origin='lower', vmin=0, vmax=8.0)
axes[0, 1].set_title("Velocidad del Reotransductor ($d\\tau/dt$)", fontweight='bold')
fig.colorbar(im_rate, ax=axes[0, 1], label="Flujo local de Presente")

# Plot 3: Tiempo Emergente Acumulado (tau)
im_tau = axes[1, 0].imshow(tau, cmap='viridis', origin='lower', vmin=0, vmax=100.0)
axes[1, 0].set_title("Tiempo Acumulado (Coordenada $\\tau$)", fontweight='bold')
fig.colorbar(im_tau, ax=axes[1, 0], label="Segundos emergentes")

# Plot 4: Integridad de los Índices de Negentropía (I)
im_index = axes[1, 1].imshow(I, cmap='cividis', origin='lower', vmin=0.0, vmax=1.0)
axes[1, 1].set_title("Integridad del Índice ($I$)", fontweight='bold')
fig.colorbar(im_index, ax=axes[1, 1], label="Orden estructural [0 a 1]")

# Recuadros y etiquetas para las estructuras A, B y C
axes[1, 1].add_patch(mpatches.Rectangle((17.5, 17.5), 4, 4, linewidth=1.5, edgecolor='#00ffff', facecolor='none'))
axes[1, 1].add_patch(mpatches.Rectangle((34.5, 34.5), 4, 4, linewidth=1.5, edgecolor='#39ff14', facecolor='none'))
axes[1, 1].add_patch(mpatches.Rectangle((39.5, 4.5), 4, 4, linewidth=1.5, edgecolor='#ff3131', facecolor='none'))

axes[1, 1].text(19.5, 22.5, 'A', color='#00ffff', fontweight='bold', fontsize=11, ha='center', va='bottom')
axes[1, 1].text(36.5, 39.5, 'B', color='#39ff14', fontweight='bold', fontsize=11, ha='center', va='bottom')
axes[1, 1].text(41.5, 9.5, 'C', color='#ff3131', fontweight='bold', fontsize=11, ha='center', va='bottom')

# Leyenda explicativa de los índices
legend_elements = [
    mpatches.Patch(facecolor='none', edgecolor='#00ffff', linewidth=1.5, label='A: Cerca Caldera 1 (1000 K) - Sostenida por alto flujo'),
    mpatches.Patch(facecolor='none', edgecolor='#39ff14', linewidth=1.5, label='B: Cerca Caldera 2 (800 K) - Sostenida por flujo medio'),
    mpatches.Patch(facecolor='none', edgecolor='#ff3131', linewidth=1.5, label='C: Aislada en vacío (2.73 K) - Decae por Landauer')
]
axes[1, 1].legend(handles=legend_elements, loc='upper left', fontsize=7.5, framealpha=0.85)

def animate(frame):
    global T, I, tau
    T, I, tau, d_tau_dt = update_physics(T, I, tau)
    
    im_temp.set_array(T)
    im_rate.set_array(d_tau_dt)
    im_tau.set_array(tau)
    im_index.set_array(I)
    
    # Ajuste suave del rango de tau conforme transcurre la simulación
    max_tau = max(1.0, np.max(tau))
    im_tau.set_clim(vmin=0, vmax=max_tau)
    
    return im_temp, im_rate, im_tau, im_index

ani = FuncAnimation(fig, animate, frames=200, interval=30, blit=False)
plt.tight_layout()
plt.show()