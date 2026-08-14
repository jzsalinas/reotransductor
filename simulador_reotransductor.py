import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button

# =====================================================================
# CONFIGURACIÓN FÍSICA INMUTABLE
# =====================================================================
GRID_SIZE = 50          
DT = 0.1                # dt FIJO: Garantiza estabilidad CFL
DIFFUSION_COEFF = 0.5   
KAPPA = 50.0            
LANDAUER_DECAY = 0.02   

# Rejillas
T = np.ones((GRID_SIZE, GRID_SIZE)) * 2.73
T[15:18, 15:18] = 1000.0
T[32:35, 32:35] = 800.0

I = np.zeros((GRID_SIZE, GRID_SIZE))
I[18:22, 18:22] = 1.0  # A
I[35:39, 35:39] = 1.0  # B
I[5:9, 40:44] = 1.0    # C

tau = np.zeros((GRID_SIZE, GRID_SIZE))
steps_per_frame = 1    # Controlado por el Slider

# =====================================================================
# MOTOR FÍSICO
# =====================================================================
def update_physics(T, I, tau):
    laplacian = np.zeros_like(T)
    laplacian[1:-1, 1:-1] = (
        T[:-2, 1:-1] + T[2:, 1:-1] +
        T[1:-1, :-2] + T[1:-1, 2:] - 4 * T[1:-1, 1:-1]
    )
    
    T_next = T + DIFFUSION_COEFF * laplacian * DT
    T_next[15:18, 15:18] = 1000.0
    T_next[32:35, 32:35] = 800.0
    T_next[0, :] = T_next[-1, :] = T_next[:, 0] = T_next[:, -1] = 2.73
    
    inv_T = 1.0 / T_next
    grad_inv_T_y, grad_inv_T_x = np.gradient(inv_T)
    
    grad_T_y, grad_T_x = np.gradient(T_next)
    J_x = -DIFFUSION_COEFF * grad_T_x
    J_y = -DIFFUSION_COEFF * grad_T_y
    
    d_tau_dt = KAPPA * np.maximum(0.0, (J_x * grad_inv_T_x + J_y * grad_inv_T_y))
    tau_next = tau + d_tau_dt * DT
    
    energy_flux = np.sqrt(J_x**2 + J_y**2)
    thermal_noise = 0.001 * T_next
    sustenance = 5.0 * energy_flux
    
    # Solo los índices vivos (I > 0) evolucionan
    mask = I > 0
    delta_I = np.zeros_like(I)
    delta_I[mask] = (sustenance[mask] - thermal_noise[mask] - LANDAUER_DECAY) * DT
    I_next = np.clip(I + delta_I, 0.0, 1.0)
    
    return T_next, I_next, tau_next, d_tau_dt

# =====================================================================
# INTERFAZ CON SLIDER
# =====================================================================
fig = plt.figure(figsize=(11, 9.5))
fig.suptitle("Simulación del Reotransductor del Presente Activo", fontsize=13, fontweight='bold')

# Subplots
ax_temp  = plt.subplot2grid((10, 2), (0, 0), rowspan=4)
ax_rate  = plt.subplot2grid((10, 2), (0, 1), rowspan=4)
ax_tau   = plt.subplot2grid((10, 2), (4, 0), rowspan=4)
ax_index = plt.subplot2grid((10, 2), (4, 1), rowspan=4)

im_temp = ax_temp.imshow(T, cmap='inferno', origin='lower', vmin=2.73, vmax=1000.0)
ax_temp.set_title("Temperatura ($T$)")
fig.colorbar(im_temp, ax=ax_temp, label="K")

im_rate = ax_rate.imshow(np.zeros((GRID_SIZE, GRID_SIZE)), cmap='plasma', origin='lower', vmin=0, vmax=8.0)
ax_rate.set_title("Velocidad Reotransductor ($d\\tau/dt$)")
fig.colorbar(im_rate, ax=ax_rate, label="Flujo Presente")

im_tau = ax_tau.imshow(tau, cmap='viridis', origin='lower', vmin=0, vmax=100.0)
ax_tau.set_title("Tiempo Acumulado ($\\tau$)")
fig.colorbar(im_tau, ax=ax_tau, label="Segundos")

im_index = ax_index.imshow(I, cmap='cividis', origin='lower', vmin=0.0, vmax=1.0)
ax_index.set_title("Integridad del Índice ($I$)")
fig.colorbar(im_index, ax=ax_index, label="Orden [0-1]")

# Espacio para el Slider
ax_slider = plt.subplot2grid((10, 2), (9, 0), colspan=2)
slider_speed = Slider(
    ax=ax_slider,
    label='Velocidad de Simulación (Pasos / Frame) ',
    valmin=1,
    valmax=50,
    valinit=1,
    valstep=1,
    color='#3b82f6'
)

def update_speed(val):
    global steps_per_frame
    steps_per_frame = int(slider_speed.val)

slider_speed.on_changed(update_speed)

# Animación con Substepping
def animate(frame):
    global T, I, tau
    
    # Ejecutamos N pasos de física según el slider antes de redibujar
    for _ in range(steps_per_frame):
        T, I, tau, d_tau_dt = update_physics(T, I, tau)
    
    im_temp.set_array(T)
    im_rate.set_array(d_tau_dt)
    im_tau.set_array(tau)
    im_index.set_array(I)
    
    max_tau = max(1.0, np.max(tau))
    im_tau.set_clim(vmin=0, vmax=max_tau)
    
    return im_temp, im_rate, im_tau, im_index

ani = FuncAnimation(fig, animate, interval=30, blit=False)
plt.tight_layout()
plt.subplots_adjust(bottom=0.08)
plt.show()