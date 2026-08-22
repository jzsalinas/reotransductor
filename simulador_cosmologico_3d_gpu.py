"""
Simulador Cosmológico 3D Interactivo Acelerado por GPU (CuPy / CUDA)
Visualizador de escritorio con Matplotlib que consume el motor central CosmologicalEngine en GPU.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
from matplotlib.colors import LinearSegmentedColormap

from server.engine import CosmologicalEngine

# Paleta oficial calibrada de la misión Planck (ESA Planck Legacy Archive)
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

# Inicializar motor central con aceleración GPU (Única Fuente de Verdad)
GRID_SIZE = 32
engine = CosmologicalEngine(grid_size=GRID_SIZE, use_gpu=True, auto_resume=True, initial_speed=20)
backend_str = "GPU (CuPy / CUDA)" if engine.gpu_available else "CPU (NumPy Fallback)"

# =====================================================================
# DASHBOARD CIENTÍFICO 3D (3x3 - 9 PANELES UNIFICADOS)
# =====================================================================
fig = plt.figure(figsize=(13.6, 9.8))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.38, wspace=0.32, top=0.93, bottom=0.08, left=0.06, right=0.94)

ax_3d      = fig.add_subplot(gs[0, 0], projection='3d')
ax_cmb     = fig.add_subplot(gs[0, 1], projection='mollweide')
ax_info    = fig.add_subplot(gs[0, 2])

ax_rho     = fig.add_subplot(gs[1, 0])
ax_rate    = fig.add_subplot(gs[1, 1])
ax_index   = fig.add_subplot(gs[1, 2])

ax_tau     = fig.add_subplot(gs[2, 0])
ax_log_tau = fig.add_subplot(gs[2, 1])
ax_temp    = fig.add_subplot(gs[2, 2])

z0 = GRID_SIZE // 2
payload = engine.get_visual_payload()
telemetry = payload["telemetry"]

# 1. Red Cósmica 3D
ax_3d.view_init(elev=24, azim=-50)
ax_3d.set_box_aspect([1, 1, 1])
sc_3d = ax_3d.scatter([], [], [], cmap='magma', s=12, alpha=0.6, vmin=0.0, vmax=4.0)
ax_3d.set_title("1. Red Cósmica 3D (Masa / Filamentos)", fontweight='bold', fontsize=8.8, pad=8)
ax_3d.set_xlim(0, GRID_SIZE)
ax_3d.set_ylim(0, GRID_SIZE)
ax_3d.set_zlim(0, GRID_SIZE)
ax_3d.tick_params(pad=0.5, labelsize=7)

# 2. CMB Mollweide
im_cmb = ax_cmb.pcolormesh(engine.LON, engine.LAT, np.array(payload["cmb"]), cmap=PLANCK_CMAP, shading='gouraud', vmin=-2.5, vmax=2.5)
ax_cmb.set_title("2. Fondo Cósmico CMB (Mollweide S² | ΔT/T̄)", fontweight='bold', fontsize=8.8, pad=8)
ax_cmb.grid(True, alpha=0.25)
cbar_cmb = fig.colorbar(im_cmb, ax=ax_cmb, orientation='horizontal', fraction=0.046, pad=0.10)
cbar_cmb.set_label("Anisotropías Relativas ΔT/T̄ (±2.5σ)", fontsize=7.2)

# 3. Consola de Telemetría
ax_info.axis('off')
ax_info.text(0.5, 0.98, f"TELEMETRÍA 3D [{backend_str}]", fontweight='bold', fontsize=9.5, ha='center', va='top', color='#0f172a')
telemetry_box = ax_info.text(
    0.5, 0.88, "", fontsize=7.4, fontfamily='monospace', va='top', ha='center',
    bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8fafc', edgecolor='#cbd5e1')
)

# 4. Densidad de Materia
im_rho = ax_rho.imshow(np.array(payload["slice_rho"]), cmap='magma', origin='lower', vmin=0.0, vmax=4.0)
title_rho = ax_rho.set_title(f"4. Densidad Materia (ρ) (z={z0})", fontweight='bold', fontsize=8.8)
div_rho = make_axes_locatable(ax_rho)
cbar_rho = fig.colorbar(im_rho, cax=div_rho.append_axes('right', size='5%', pad=0.06))
cbar_rho.set_label("Densidad [ρ/ρ̄]", fontsize=7.2)

# 5. Velocidad Reotransductor
im_rate = ax_rate.imshow(np.array(payload["slice_rate"]), cmap='plasma', origin='lower', vmin=0.0, vmax=1.5)
title_rate = ax_rate.set_title(f"5. Velocidad Reotransductor (dτ/dt) (z={z0})", fontweight='bold', fontsize=8.8)
div_rate = make_axes_locatable(ax_rate)
cbar_rate = fig.colorbar(im_rate, cax=div_rate.append_axes('right', size='5%', pad=0.06))
cbar_rate.set_label("Flujo Local de Presente", fontsize=7.2)

# 6. Autoorganización
im_index = ax_index.imshow(np.array(payload["slice_index"]), cmap='cividis', origin='lower', vmin=0.0, vmax=1.0)
title_index = ax_index.set_title(f"6. Autoorganización I(r,t) (z={z0})", fontweight='bold', fontsize=8.8)
div_index = make_axes_locatable(ax_index)
cbar_index = fig.colorbar(im_index, cax=div_index.append_axes('right', size='5%', pad=0.06))
cbar_index.set_label("Orden Informacional [0 a 1]", fontsize=7.2)

# 7. Tiempo Lineal
im_tau = ax_tau.imshow(np.array(payload["slice_tau"]), cmap='viridis', origin='lower', vmin=0.0, vmax=100.0)
title_tau = ax_tau.set_title(f"7. Tiempo Lineal (τ) (z={z0})", fontweight='bold', fontsize=8.8)
div_tau = make_axes_locatable(ax_tau)
cbar_tau = fig.colorbar(im_tau, cax=div_tau.append_axes('right', size='5%', pad=0.06))
cbar_tau.set_label("Segundos Propios τ", fontsize=7.2)

# 8. Red Fósil
im_log_tau = ax_log_tau.imshow(np.array(payload["slice_log_tau"]), cmap='inferno', origin='lower', vmin=0.0, vmax=4.0)
title_log_tau = ax_log_tau.set_title(f"8. Red Fósil log₁₀(1+τ) (z={z0})", fontweight='bold', fontsize=8.8)
div_log_tau = make_axes_locatable(ax_log_tau)
cbar_log_tau = fig.colorbar(im_log_tau, cax=div_log_tau.append_axes('right', size='5%', pad=0.06))
cbar_log_tau.set_label("Memoria Arqueológica", fontsize=7.2)

# 9. Temperatura
im_temp = ax_temp.imshow(np.array(payload["slice_temp"]), cmap='inferno', origin='lower', vmin=2.73, vmax=45.0)
title_temp = ax_temp.set_title(f"9. Temperatura Plasma T (z={z0})", fontweight='bold', fontsize=8.8)
div_temp = make_axes_locatable(ax_temp)
cbar_temp = fig.colorbar(im_temp, cax=div_temp.append_axes('right', size='5%', pad=0.06))
cbar_temp.set_label("Kelvins [K] (0°C = 273.15 K)", fontsize=7.2)

# Slider de Velocidad (hasta 1000x)
ax_slider = fig.add_axes([0.22, 0.022, 0.56, 0.022])
slider_speed = Slider(
    ax=ax_slider,
    label='Aceleración (Pasos / Frame) ',
    valmin=1,
    valmax=1000,
    valinit=engine.steps_per_frame,
    valstep=1,
    color='#2563eb'
)

def update_speed(val):
    engine.steps_per_frame = int(slider_speed.val)

slider_speed.on_changed(update_speed)


def animate_3d(frame):
    # 1. Integración en GPU en el motor central
    engine.step_batch()

    # 2. Obtener datos serializados
    payload = engine.get_visual_payload()
    tel = payload["telemetry"]
    z_slice = tel["z_slice"]

    # 3. Puntos 3D
    pts = payload["points_3d"]
    if len(pts) > 0:
        pts_arr = np.array(pts)
        sc_3d._offsets3d = (pts_arr[:, 0], pts_arr[:, 1], pts_arr[:, 2])
        sc_3d.set_array(pts_arr[:, 3])
        sc_3d.set_clim(vmin=0.0, vmax=max(3.0, float(np.max(pts_arr[:, 3]))))

    # 4. CMB
    cmb_data = np.array(payload["cmb"])
    im_cmb.set_array(cmb_data.ravel())

    # 5. Cortes 2D
    im_rho.set_array(np.array(payload["slice_rho"]))
    im_rate.set_array(np.array(payload["slice_rate"]))
    im_index.set_array(np.array(payload["slice_index"]))
    im_tau.set_array(np.array(payload["slice_tau"]))
    im_log_tau.set_array(np.array(payload["slice_log_tau"]))
    im_temp.set_array(np.array(payload["slice_temp"]))

    # Títulos y Telemetría
    fig.suptitle(f"Simulación Cosmológica 3D [GPU]: Reotransductor, CMB & Bekenstein | Eón N = {tel['eon']} [{tel['era']}]", fontsize=11.5, fontweight='bold')
    title_rho.set_text(f"4. Densidad Materia (ρ) [a={tel['scale_factor']:.2f}] (Plano Z = {z_slice})")
    title_rate.set_text(f"5. Velocidad Reotransductor (dτ/dt) (Plano Z = {z_slice})")
    title_index.set_text(f"6. Autoorganización I(r,t) (Plano Z = {z_slice})")
    title_tau.set_text(f"7. Tiempo Lineal (τ) (Plano Z = {z_slice})")
    title_log_tau.set_text(f"8. Red Fósil log₁₀(1+τ) (Plano Z = {z_slice})")
    title_temp.set_text(f"9. Temperatura Plasma T (Plano Z = {z_slice})")

    bar_len = 10
    filled_len = int((tel["tunnel_progress"] / 100.0) * bar_len)
    bar_str = "█" * filled_len + "░" * (bar_len - filled_len)

    telemetry_str = (
        f"• Motor de Cómputo: {backend_str}\n"
        f"• Eón Cósmico Actual: N = {tel['eon']}\n"
        f"• Factor de Escala: a = {tel['scale_factor']:.3f} (z = {tel['redshift']:.2f})\n"
        f"• Temperatura Plasma: {tel['temp_norm']:.1f} K ({tel['temp_astro']:.0f} K Astrofísico)\n"
        f"• Velocidad de la Luz (c): {tel['c_light']:.2f} celdas/s\n"
        f"• Masa Núcleo (M_BH): {tel['mass_fraction']:.1f}% del total\n"
        f"• Atractor Central 3D: (X={tel['attractor']['x']}, Y={tel['attractor']['y']}, Z={tel['attractor']['z']})\n"
        f"• Entropía Eón N (S_BH): {tel['s_bh']:.0f} k_B\n"
        f"• Límite Bekenstein (S_max): {tel['s_crit']:.0f} k_B\n"
        f"• {tel['progress_label']}: [{bar_str}] {tel['tunnel_progress']:.1f}%\n"
        f"• Odómetro Fósil Total 3D: {tel['fossil_odometer']:.0f} s ({tel['time_myr']:.1f} Myr / {tel['time_myr']/1000.0:.2f} Gyr)\n"
        f"• Estado: {tel['state_status']}"
    )
    telemetry_box.set_text(telemetry_str)

    return sc_3d, im_cmb, im_rho, im_rate, im_index, im_tau, im_log_tau, im_temp, telemetry_box


if __name__ == '__main__':
    ani = FuncAnimation(fig, animate_3d, interval=30, blit=False, cache_frame_data=False)
    plt.show()
