import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# =====================================================================
# CONFIGURACIÓN DE PARÁMETROS FÍSICOS Y GEOMETRÍA
# =====================================================================
GRID_SIZE = 50          # Resolución de la Malla espacial
DT = 0.1                # Intervalo de integración temporal (dt)
DIFFUSION_COEFF = 0.5   # Conductividad térmica del medio (k)
KAPPA = 50.0            # Constante disipativa del Reotransductor (kappa)
LANDAUER_DECAY = 0.02   # Tasa de decaimiento entrópico natural de los Índices

# Inicialización de rejillas físicas
# 1. Temperatura (T): Inicializada a la temperatura del vacío frío (2.73 K)
T = np.ones((GRID_SIZE, GRID_SIZE)) * 2.73

# Colocar fuentes termodinámicas calientes (Estrellas / Calderas a 1000 K)
T[15:18, 15:18] = 1000.0
T[32:35, 32:35] = 800.0

# 2. Índices de Baja Entropía (I): Estructuras ordenadas iniciales (0 a 1)
# Colocamos "islas de orden" (como sistemas vivos o bases de datos) en la malla
I = np.zeros((GRID_SIZE, GRID_SIZE))
I[18:22, 18:22] = 1.0  # Estructura A (Cerca de la caldera 1)
I[35:39, 35:39] = 1.0  # Estructura B (Cerca de la caldera 2)
I[5:9, 40:44] = 1.0    # Estructura C (Aislada en el vacío absoluto)

# 3. Tiempo Térmico Acumulado (tau): El reloj local emergente de cada celda
tau = np.zeros((GRID_SIZE, GRID_SIZE))

# =====================================================================
# MOTOR DE EVOLUCIÓN FÍSICA
# =====================================================================
def update_physics(T, I, tau):
    # a. Difusión de Calor (Ecuación de conducción de Fourier)
    # Laplaciano de Temperatura usando diferencias finitas
    laplacian = (
        np.roll(T, 1, axis=0) + np.roll(T, -1, axis=0) +
        np.roll(T, 1, axis=1) + np.roll(T, -1, axis=1) - 4 * T
    )
    # Mantener las calderas de energía constantes (focos calientes persistentes)
    T_next = T + DIFFUSION_COEFF * laplacian * DT
    T_next[15:18, 15:18] = 1000.0
    T_next[32:35, 32:35] = 800.0
    # Fronteras frías (Disipación hacia el espacio profundo a 2.73 K)
    T_next[0, :] = T_next[-1, :] = T_next[:, 0] = T_next[:, -1] = 2.73
    
    # b. Gradiente de Temperatura Inversa (Fuerza impulsora: X = grad(1/T))
    inv_T = 1.0 / T_next
    grad_inv_T_y, grad_inv_T_x = np.gradient(inv_T)
    
    # c. Corriente disipativa (Flujo de energía J = -k * grad(T))
    grad_T_y, grad_T_x = np.gradient(T_next)
    J_x = -DIFFUSION_COEFF * grad_T_x
    J_y = -DIFFUSION_COEFF * grad_T_y
    
    # d. Velocidad del Reotransductor (d_tau/dt = kappa * |J . grad(1/T)|)
    # Es el producto escalar local de la disipación de Onsager
    d_tau_dt = KAPPA * np.abs(J_x * grad_inv_T_x + J_y * grad_inv_T_y)
    
    # Actualizar la coordenada del "Presente Activo" (tau)
    tau_next = tau + d_tau_dt * DT
    
    # e. Dinámica de los Índices (I) bajo el Límite de Landauer
    # Para sobrevivir, un Índice necesita consumir un flujo de Joules (J) proporcional
    # a la disipación local para compensar el desorden térmico de fondo (T)
    energy_flux = np.sqrt(J_x**2 + J_y**2)
    
    I_next = np.copy(I)
    for y in range(GRID_SIZE):
        for x in range(GRID_SIZE):
            if I[y, x] > 0:
                # El calor excesivo (ruido térmico) desorganiza el índice
                thermal_noise = 0.001 * T_next[y, x]
                # El flujo local de energía disipada sostiene y alimenta al índice (negentropía)
                sustenance = 5.0 * energy_flux[y, x]
                
                # Ecuación de balance para la integridad del Índice
                delta_I = (sustenance - thermal_noise - LANDAUER_DECAY) * DT
                I_next[y, x] = np.clip(I[y, x] + delta_I, 0.0, 1.0)
                
    return T_next, I_next, tau_next, d_tau_dt

# =====================================================================
# INTERFAZ GRÁFICA DE VISUALIZACIÓN
# =====================================================================
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
fig.suptitle("Simulación del Reotransductor del Presente Activo", fontsize=15, fontweight='bold')

# Plot 1: Temperatura
im_temp = axes[0, 0].imshow(T, cmap='inferno', origin='lower')
axes[0, 0].set_title("Temperatura de la Malla ($T$)", fontsize=11, fontweight='bold')
fig.colorbar(im_temp, ax=axes[0, 0], label="Kelvins [K]")

# Plot 2: Velocidad de flujo del tiempo (d_tau/dt)
im_rate = axes[0, 1].imshow(np.zeros((GRID_SIZE, GRID_SIZE)), cmap='plasma', origin='lower')
axes[0, 1].set_title("Velocidad del Reotransductor ($d\\tau/dt$)", fontsize=11, fontweight='bold')
fig.colorbar(im_rate, ax=axes[0, 1], label="Flujo local de Presente")

# Plot 3: Tiempo Emergente Acumulado (tau)
im_tau = axes[1, 0].imshow(tau, cmap='viridis', origin='lower')
axes[1, 0].set_title("Tiempo Acumulado (Coordenada $\\tau$)", fontsize=11, fontweight='bold')
fig.colorbar(im_tau, ax=axes[1, 0], label="Segundos emergentes")

# Plot 4: Integridad de los Índices de Negentropía (I)
im_index = axes[1, 1].imshow(I, cmap='Blues_r', origin='lower')
axes[1, 1].set_title("Integridad del Índice ($I$)", fontsize=11, fontweight='bold')
fig.colorbar(im_index, ax=axes[1, 1], label="Estabilidad estructural")

def animate(frame):
    global T, I, tau
    T, I, tau, d_tau_dt = update_physics(T, I, tau)
    
    im_temp.set_array(T)
    im_rate.set_array(d_tau_dt)
    im_tau.set_array(tau)
    im_index.set_array(I)
    
    # Ajustar escalas dinámicamente para la tasa de tiempo
    im_rate.set_clim(vmin=0, vmax=max(1e-3, np.max(d_tau_dt)))
    im_tau.set_clim(vmin=0, vmax=max(1e-3, np.max(tau)))
    
    return im_temp, im_rate, im_tau, im_index

ani = FuncAnimation(fig, animate, frames=200, interval=50, blit=False)
plt.tight_layout()
plt.show()