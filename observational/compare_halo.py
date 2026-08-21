"""
Command-Line Utility and Publication Plot Generator:
Compares Reotransductor Cosmological Simulation with SPARC Galaxy Observations and NFW Cusp Model.
Evaluates the resolution of the Cusp-Core Problem and generates assets/halo_cusp_core_sparc.png.
"""

import os
import glob
import shutil
import argparse
from typing import Dict, Any, Optional
import matplotlib.pyplot as plt
import numpy as np

from server.engine import CosmologicalEngine
from observational.halo_data import SPARCHaloData
from observational.halo_analyzer import HaloRadialProfileAnalyzer


def run_halo_comparison(
    checkpoint_path: str = "checkpoints",
    auto_resume: bool = True,
    steps: int = 0,
    galaxy_name: str = "DDO_154",
    output_fig: str = "assets/halo_cusp_core_sparc.png"
) -> Dict[str, Any]:
    """Executes live halo radial profile comparison on active cosmological state."""
    os.makedirs(os.path.dirname(output_fig) or ".", exist_ok=True)
    print("=" * 75)
    print("  🌌 REOTRANSDUCTOR OBSERVATIONAL PIPELINE: CUSP-CORE & SPARC VALIDATION")
    print("=" * 75)
    print(f"• Loading Cosmological State (checkpoint_dir='{checkpoint_path}', auto_resume={auto_resume})...")
    engine = CosmologicalEngine(checkpoint_dir=checkpoint_path, auto_resume=auto_resume)
    if steps > 0:
        print(f"• Evolving additional {steps} integration steps...")
        for _ in range(steps):
            engine.step()

    print(f"• Evaluated State: Eon {engine.eon} | Steps: {engine.total_steps:,} | Scale Factor a = {engine.scale_factor:.3f}")
    metrics = generate_eon_halo_report(engine, galaxy_name=galaxy_name, output_dir=os.path.dirname(output_fig) or ".")
    print("\n" + "=" * 75)
    print("  📊 DARK MATTER HALO & CUSP-CORE DIAGNOSTICS")
    print("=" * 75)
    print(f"• Halo Attractor Center:          ({metrics['halo_center']['x']}, {metrics['halo_center']['y']}, {metrics['halo_center']['z']})")
    print(f"• Inner Logarithmic Slope gamma_0: {metrics['inner_slope_gamma0']:.3f} ({'✅ FLAT CORE' if metrics['is_cored'] else '⚠️ CUSPY'})")
    print(f"• Model Preference:               {metrics['fit_results']['preferred_model']}")
    print(f"• Chi^2 Burkert Core vs NFW Cusp: {metrics['fit_results']['chi2_core']} (Core) vs {metrics['fit_results']['chi2_nfw']} (NFW)")
    print(f"• Output Plot:                     {output_fig}")
    print("=" * 75)
    return metrics


def generate_eon_halo_report(
    engine,
    galaxy_name: str = "DDO_154",
    output_dir: str = "checkpoints/snapshots"
) -> Dict[str, Any]:
    """
    Generates complete Halo & Cusp-Core observational report for active or checkpointed eon:
      1. Spherically averages 3D density profile around primary halo center
      2. Computes logarithmic slope gamma(r) = d ln(rho) / d ln(r)
      3. Computes rotation curve V_c(r)
      4. Renders and saves 3-panel publication figure
    """
    os.makedirs(output_dir, exist_ok=True)
    eon = engine.eon
    scale_factor = float(engine.scale_factor)

    # 1. Ingest SPARC Galaxy Observational Data
    sparc = SPARCHaloData()
    gal_data = sparc.get_galaxy(galaxy_name)

    # 2. Extract Halo Profiles from 3D Simulation
    analyzer = HaloRadialProfileAnalyzer(
        grid_size=engine.grid_size,
        box_size_mpc=100.0,
        n_shells=24
    )
    
    # Gravitational potential if available
    phi_3d = getattr(engine, 'phi', None)
    halo_metrics = analyzer.evaluate_halo_diagnostics(
        rho_3d=engine.rho,
        phi_3d=phi_3d,
        tau_3d=engine.tau
    )

    r_sim = np.array(halo_metrics["r_mpc"])
    rho_sim = np.array(halo_metrics["rho_radial"])
    gamma_sim = np.array(halo_metrics["log_slope_gamma"])
    v_sim = np.array(halo_metrics["v_circular"])

    # Physical galactic rescaling (normalize to SPARC galaxy radius and velocity scale)
    r_max_gal = float(np.max(gal_data["r_kpc"]))
    v_max_gal = float(gal_data["v_flat_kms"])
    
    # Scale factors for physical overlay
    scale_r = r_max_gal / np.max(r_sim)
    r_sim_kpc = r_sim * scale_r
    
    v_sim_max = max(1e-3, np.max(v_sim))
    v_sim_kms = (v_sim / v_sim_max) * v_max_gal

    rho_sim_max = max(1e-3, np.max(rho_sim))
    rho_gal_max = max(1e-4, np.max(gal_data["rho_dm"])) if len(gal_data["rho_dm"]) > 0 else 0.04
    rho_sim_norm = (rho_sim / rho_sim_max) * rho_gal_max

    # Theoretical Benchmark Models
    r_theory = np.linspace(0.1, r_max_gal, 100)
    r_s_fit = float(halo_metrics["fit_results"].get("fitted_scale_radius_nfw_mpc", 8.0)) * scale_r
    r_c_fit = float(halo_metrics["fit_results"].get("fitted_core_radius_mpc", 6.0)) * scale_r

    rho_nfw_theory = sparc.nfw_density(r_theory, rho_s=rho_gal_max * 0.4, r_s=r_s_fit)
    rho_core_theory = sparc.burkert_density(r_theory, rho_0=rho_gal_max, r_0=r_c_fit)

    gamma_nfw_theory = sparc.nfw_log_slope(r_theory, r_s=r_s_fit)
    gamma_core_theory = sparc.burkert_log_slope(r_theory, r_0=r_c_fit)

    # 3. Render 3-Panel Publication Figure
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(10, 11),
        gridspec_kw={'height_ratios': [1.8, 1.2, 1.4], 'hspace': 0.28}
    )
    fig.patch.set_facecolor('#0d131f')

    for ax in (ax1, ax2, ax3):
        ax.set_facecolor('#131b2e')
        ax.grid(True, linestyle=':', alpha=0.35, color='#475569')
        ax.tick_params(colors='#94a3b8', labelsize=10)
        for spine in ax.spines.values():
            spine.set_edgecolor('#334155')

    # -------------------------------------------------------------------------
    # Panel 1: Radial Density Profile rho(r)
    # -------------------------------------------------------------------------
    # NFW Cusp Theory
    ax1.plot(
        r_theory, rho_nfw_theory,
        color='#ef4444', linestyle='--', linewidth=2.0,
        label=r'Standard $\Lambda$CDM NFW Cusp ($\rho \propto r^{-1}$ as $r \to 0$)'
    )

    # Burkert / Reotransductor Core Theory
    ax1.plot(
        r_theory, rho_core_theory,
        color='#06b6d4', linestyle='-.', linewidth=2.0,
        label=r'Burkert Cored Model ($\rho \to \mathrm{const}$ as $r \to 0$)'
    )

    # SPARC Galaxy Data Points if density available
    if len(gal_data["rho_dm"]) > 0 and np.any(gal_data["rho_dm"] > 0):
        ax1.scatter(
            gal_data["r_kpc"], gal_data["rho_dm"],
            color='#38bdf8', marker='o', s=40, edgecolors='#0284c7', zorder=5,
            label=f'SPARC {gal_data["name"]} Inferred DM Density'
        )

    # Reotransductor Simulation Profile
    ax1.plot(
        r_sim_kpc, rho_sim_norm,
        color='#f59e0b', marker='s', markersize=5, linewidth=2.4, zorder=6,
        label=f'Reotransductor Eon {eon} Virialized Core (3D Hydro)'
    )

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(0.2, r_max_gal * 1.1)
    ax1.set_ylabel(r'Density $\rho(r)\ \ [M_\odot/\mathrm{pc}^3]$', color='#f8fafc', fontsize=11, fontweight='bold')
    ax1.set_title(
        f'Dark Matter Halo Profile & Cusp-Core Problem — Eon {eon} vs. SPARC {gal_data["name"]}\n'
        f'(Inner Logarithmic Slope $\\gamma_0 = {halo_metrics["inner_slope_gamma0"]:.3f}$ | Preferred: {halo_metrics["fit_results"]["preferred_model"]})',
        color='#f8fafc', fontsize=12, fontweight='bold', pad=10
    )
    leg1 = ax1.legend(loc='upper right', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=9)
    for t in leg1.get_texts():
        t.set_color('#e2e8f0')

    # -------------------------------------------------------------------------
    # Panel 2: Logarithmic Density Slope gamma(r) = d ln(rho) / d ln(r)
    # -------------------------------------------------------------------------
    ax2.axhline(0.0, color='#06b6d4', linestyle=':', linewidth=1.5, label='Flat Core Boundary (gamma = 0.0)')
    ax2.axhline(-1.0, color='#ef4444', linestyle=':', linewidth=1.5, label='NFW Cusp Boundary (gamma = -1.0)')
    ax2.axhline(-3.0, color='#94a3b8', linestyle=':', linewidth=1.2, label='Outer Envelope (gamma = -3.0)')

    ax2.plot(r_theory, gamma_nfw_theory, color='#ef4444', linestyle='--', linewidth=1.8)
    ax2.plot(r_theory, gamma_core_theory, color='#06b6d4', linestyle='-.', linewidth=1.8)

    ax2.plot(
        r_sim_kpc, gamma_sim,
        color='#f59e0b', marker='o', markersize=4.5, linewidth=2.2,
        label='Simulation Slope $\\gamma(r) = \\mathrm{d}\\ln\\rho / \\mathrm{d}\\ln r$'
    )

    ax2.set_xscale('log')
    ax2.set_xlim(0.2, r_max_gal * 1.1)
    ax2.set_ylim(-3.5, 0.5)
    ax2.set_ylabel(r'Log Slope $\gamma(r)$', color='#f8fafc', fontsize=10, fontweight='bold')
    leg2 = ax2.legend(loc='lower left', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=8.5)
    for t in leg2.get_texts():
        t.set_color('#e2e8f0')

    # -------------------------------------------------------------------------
    # Panel 3: Circular Velocity Rotation Curve V_c(r)
    # -------------------------------------------------------------------------
    # SPARC Observed Rotation Curve with Error Bars
    ax3.errorbar(
        gal_data["r_kpc"], gal_data["v_obs_kms"], yerr=gal_data["err_v_kms"],
        fmt='o', color='#38bdf8', ecolor='#0284c7', elinewidth=1.6, capsize=3.5,
        markersize=6, label=f'SPARC {gal_data["name"]} Observed $V_{{\\mathrm{{obs}}}}(r)$'
    )

    # Reotransductor Circular Velocity
    ax3.plot(
        r_sim_kpc, v_sim_kms,
        color='#f59e0b', marker='s', markersize=4.5, linewidth=2.2,
        label=f'Reotransductor $V_c(r) = \\sqrt{{GM(<r)/r}}$'
    )

    ax3.set_xlim(0.0, r_max_gal * 1.05)
    ax3.set_ylim(0.0, v_max_gal * 1.35)
    ax3.set_xlabel('Galactocentric Radius $r\\ \\ [\\mathrm{kpc}]$', color='#f8fafc', fontsize=11, fontweight='bold')
    ax3.set_ylabel(r'Velocity $V_c\ \ [\mathrm{km/s}]$', color='#f8fafc', fontsize=10, fontweight='bold')
    leg3 = ax3.legend(loc='lower right', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=9)
    for t in leg3.get_texts():
        t.set_color('#e2e8f0')

    # Save Eon Plot
    eon_halo_path = os.path.join(output_dir, f"halo_comparison_eon_{eon}.png")
    plt.savefig(eon_halo_path, dpi=180, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)

    # Copy to assets/
    latest_halo_path = "assets/halo_cusp_core_sparc.png"
    try:
        os.makedirs("assets", exist_ok=True)
        shutil.copyfile(eon_halo_path, latest_halo_path)
    except Exception:
        pass

    return {
        "eon": eon,
        "scale_factor": scale_factor,
        "halo_center": halo_metrics["halo_center"],
        "inner_slope_gamma0": halo_metrics["inner_slope_gamma0"],
        "is_cored": halo_metrics["is_cored"],
        "fit_results": halo_metrics["fit_results"],
        "halo_comparison_path": eon_halo_path
    }


def process_all_existing_checkpoints_halo(checkpoints_dir: str = "checkpoints", galaxy_name: str = "DDO_154"):
    """Processes all historical cluster/virialized checkpoints (clusters_eon_*.npz or eon_*.npz)."""
    # Prioritize dedicated virialized cluster epoch checkpoints (a ~ 3.0)
    npz_files = sorted(glob.glob(os.path.join(checkpoints_dir, "clusters_eon_*.npz")))
    if not npz_files:
        npz_files = sorted(glob.glob(os.path.join(checkpoints_dir, "eon_*.npz")))

    if not npz_files:
        print(f"No cluster or eon checkpoints found in '{checkpoints_dir}'.")
        return

    print(f"• Found {len(npz_files)} Halo/Cluster checkpoints in '{checkpoints_dir}'. Processing Cusp-Core reports...")
    snapshots_dir = os.path.join(checkpoints_dir, "snapshots")
    engine = CosmologicalEngine(checkpoint_dir=checkpoints_dir, auto_resume=False)

    for npz_path in npz_files:
        try:
            data = np.load(npz_path)
            engine.rho = data['rho']
            engine.v_x = data['v_x']
            engine.v_y = data['v_y']
            engine.v_z = data['v_z']
            engine.T = data['T']
            engine.I = data['I']
            engine.tau = data['tau']
            if 'phi' in data:
                engine.phi = data['phi']
            engine.scale_factor = float(data['scale_factor'])
            engine.eon = int(data['eon'])
            engine.total_steps = int(data.get('total_steps', 0))
            if 'tau_eon_start' in data:
                engine.tau_eon_start = data['tau_eon_start']

            print(f"\n--- Processing Halo: {os.path.basename(npz_path)} (Eon {engine.eon}) ---")
            metrics = generate_eon_halo_report(engine, galaxy_name=galaxy_name, output_dir=snapshots_dir)
            print(f"  * Slope gamma_0: {metrics['inner_slope_gamma0']} | Core: {metrics['is_cored']} | Pref: {metrics['fit_results']['preferred_model']}")
        except Exception as ex:
            print(f"  * Error processing {npz_path}: {ex}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Reotransductor vs. SPARC Cusp-Core Dark Matter Halo Validation")
    parser.add_argument("--checkpoint-dir", type=str, default="checkpoints", help="Directory with checkpoints")
    parser.add_argument("--steps", type=int, default=0, help="Additional steps to simulate before evaluation")
    parser.add_argument("--from-scratch", action="store_true", help="Start from blank initial conditions")
    parser.add_argument("--process-all", action="store_true", help="Process and generate Halo plots for all historical checkpoints")
    parser.add_argument("--galaxy", type=str, default="DDO_154", help="SPARC benchmark galaxy identifier (e.g. DDO_154, NGC_2403)")
    parser.add_argument("--output", type=str, default="assets/halo_cusp_core_sparc.png", help="Output PNG path")
    args = parser.parse_args()

    if args.process_all:
        process_all_existing_checkpoints_halo(checkpoints_dir=args.checkpoint_dir, galaxy_name=args.galaxy)
    else:
        run_halo_comparison(
            checkpoint_path=args.checkpoint_dir,
            auto_resume=not args.from_scratch,
            steps=args.steps,
            galaxy_name=args.galaxy,
            output_fig=args.output
        )
