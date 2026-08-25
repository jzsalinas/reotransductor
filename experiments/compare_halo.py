"""
Command-Line Utility and Publication Plot Generator:
Compares Reotransductor Cosmological Simulation with SPARC Galaxy Observations and NFW Cusp Model.
Evaluates the resolution of the Cusp-Core Problem and generates results/paper_II/halo_cusp_core_sparc.png.
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
    galaxy_name: str = "NGC_2403",
    output_fig: str = "results/paper_II/halo_cusp_core_sparc.png"
) -> Dict[str, Any]:
    """Executes live halo radial profile comparison on active cosmological state."""
    os.makedirs(os.path.dirname(output_fig) or ".", exist_ok=True)
    print("=" * 75)
    print("  🌌 REOTRANSDUCTOR OBSERVATIONAL PIPELINE: CUSP-CORE & SPARC VALIDATION")
    print("=" * 75)
    
    # If checkpoint_path is a directory, prioritize virialized cluster checkpoints (a = 3.0)
    target_checkpoint = None
    if os.path.isdir(checkpoint_path):
        cluster_files = sorted(glob.glob(os.path.join(checkpoint_path, "clusters_eon_*.npz")))
        if cluster_files:
            target_checkpoint = cluster_files[-1]
    elif os.path.isfile(checkpoint_path):
        target_checkpoint = checkpoint_path

    if target_checkpoint and os.path.exists(target_checkpoint):
        print(f"• Loading Virialized Cluster Epoch Checkpoint: '{target_checkpoint}'...")
        data = np.load(target_checkpoint)
        grid_n = int(data['rho'].shape[0])
        engine = CosmologicalEngine.from_checkpoint(target_checkpoint)
    else:
        print(f"• Loading Cosmological State (checkpoint_dir='{checkpoint_path}', auto_resume={auto_resume})...")
        engine = CosmologicalEngine(checkpoint_dir=checkpoint_path, auto_resume=auto_resume)

    if steps > 0:
        print(f"• Evolving additional {steps} integration steps...")
        for _ in range(steps):
            engine.step()

    print(f"• Evaluated State: Eon {engine.eon} | Steps: {engine.total_steps:,} | Scale Factor a = {engine.scale_factor:.3f}")
    snapshots_dir = "checkpoints/snapshots"
    metrics = generate_eon_halo_report(engine, galaxy_name=galaxy_name, output_dir=snapshots_dir)
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
    galaxy_name: str = "NGC_2403",
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

    # 2. Extract Halo Radial Diagnostics from Simulation
    box_sz = getattr(engine, 'box_size_mpc', 100.0)
    analyzer = HaloRadialProfileAnalyzer(
        grid_size=engine.grid_size,
        box_size_mpc=box_sz,
        n_shells=24,
        r_max_mpc=min(90.0, box_sz / 2.0)
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

    # Copy to results/paper_II/
    latest_halo_path = "results/paper_II/halo_cusp_core_sparc.png"
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


def run_full_sparc_population_analysis(
    checkpoint_path: str = "checkpoints",
    output_fig: str = "results/paper_II/halo_cusp_core_sparc.png"
) -> Dict[str, Any]:
    """
    Executes full population-level benchmark across all 175 SPARC galaxies:
      1. Evaluates galaxy-by-galaxy fits for Burkert Core vs. NFW Cusp
      2. Computes empirical stacked rotation curve (median and 1sigma dispersion)
      3. Overlays Reotransductor emergent 3D hydrodynamic halo rotation curve
      4. Renders a comprehensive 4-panel publication figure
    """
    os.makedirs(os.path.dirname(output_fig) or ".", exist_ok=True)
    print("=" * 75)
    print("  🌌 REOTRANSDUCTOR OBSERVATIONAL PIPELINE: FULL SPARC (175 GALAXIES) POPULATION")
    print("=" * 75)
    
    # 1. Ingest and analyze all 175 SPARC galaxies
    sparc = SPARCHaloData()
    pop_stats = sparc.evaluate_full_catalog()
    
    print(f"• Total SPARC Galaxies Evaluated: {pop_stats['total_galaxies_evaluated']}")
    print(f"• Cored Profile Preference (Reotransductor): {pop_stats['core_preferred_count']}/{pop_stats['total_galaxies_evaluated']} ({pop_stats['core_preference_pct']}%)")
    print(f"• Cuspy Profile Preference (NFW):           {pop_stats['nfw_preferred_count']}/{pop_stats['total_galaxies_evaluated']} ({pop_stats['nfw_preference_pct']}%)")
    print(f"• Mean Reduced Chi^2:                       Core = {pop_stats['mean_reduced_chi2_core']} | NFW = {pop_stats['mean_reduced_chi2_nfw']}")

    # 2. Extract active simulation halo from cluster epoch
    target_checkpoint = None
    if os.path.isdir(checkpoint_path):
        cluster_files = sorted(glob.glob(os.path.join(checkpoint_path, "clusters_eon_*.npz")))
        if cluster_files:
            target_checkpoint = cluster_files[-1]
    elif os.path.isfile(checkpoint_path):
        target_checkpoint = checkpoint_path

    if target_checkpoint and os.path.exists(target_checkpoint):
        data = np.load(target_checkpoint)
        grid_n = int(data['rho'].shape[0])
        engine = CosmologicalEngine.from_checkpoint(target_checkpoint)
    else:
        engine = CosmologicalEngine(checkpoint_dir=checkpoint_path, auto_resume=True)

    box_sz = getattr(engine, 'box_size_mpc', 100.0)
    analyzer = HaloRadialProfileAnalyzer(grid_size=engine.grid_size, box_size_mpc=box_sz, n_shells=24, r_max_mpc=min(90.0, box_sz / 2.0))
    halo_metrics = analyzer.evaluate_halo_diagnostics(rho_3d=engine.rho, phi_3d=getattr(engine, 'phi', None), tau_3d=engine.tau)

    # 3. Reference Galaxy NGC 2403 Data
    gal_data = sparc.get_galaxy("NGC_2403")
    r_max_gal = float(np.max(gal_data["r_kpc"]))
    v_max_gal = float(gal_data["v_flat_kms"])
    rho_gal_max = float(gal_data.get("best_fit_central_density_msun_pc3", 0.045))

    r_sim = np.array(halo_metrics["r_mpc"])
    rho_sim = np.array(halo_metrics["rho_radial"])
    gamma_sim = np.array(halo_metrics["log_slope_gamma"])
    v_sim = np.array(halo_metrics["v_circular"])

    scale_r = r_max_gal / max(1e-3, np.max(r_sim))
    r_sim_kpc = r_sim * scale_r
    v_sim_kms = (v_sim / max(1e-3, np.max(v_sim))) * v_max_gal
    rho_sim_norm = (rho_sim / max(1e-3, np.max(rho_sim))) * rho_gal_max

    r_theory = np.linspace(0.1, r_max_gal, 100)
    r_s_fit = float(halo_metrics["fit_results"].get("fitted_scale_radius_nfw_mpc", 8.0)) * scale_r
    r_c_fit = float(halo_metrics["fit_results"].get("fitted_core_radius_mpc", 6.0)) * scale_r
    rho_nfw_theory = sparc.nfw_density(r_theory, rho_s=rho_gal_max * 0.4, r_s=r_s_fit)
    rho_core_theory = sparc.burkert_density(r_theory, rho_0=rho_gal_max, r_0=r_c_fit)
    gamma_nfw_theory = sparc.nfw_log_slope(r_theory, r_s=r_s_fit)
    gamma_core_theory = sparc.burkert_log_slope(r_theory, r_0=r_c_fit)

    # 4. Render 4-Panel Publication Figure
    fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(
        2, 2, figsize=(14, 11),
        gridspec_kw={'hspace': 0.28, 'wspace': 0.22}
    )
    fig.patch.set_facecolor('#0d131f')

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_facecolor('#131b2e')
        ax.grid(True, linestyle=':', alpha=0.35, color='#475569')
        ax.tick_params(colors='#94a3b8', labelsize=10)
        for spine in ax.spines.values():
            spine.set_edgecolor('#334155')

    # Panel 1: Radial Density Profile (NGC 2403)
    ax1.plot(r_theory, rho_nfw_theory, color='#ef4444', linestyle='--', linewidth=2.0, label=r'NFW Cusp ($\rho \propto r^{-1}$)')
    ax1.plot(r_theory, rho_core_theory, color='#06b6d4', linestyle='-.', linewidth=2.0, label=r'Burkert Core ($\rho \to \mathrm{const}$)')
    ax1.scatter(gal_data["r_kpc"], gal_data["rho_dm"], color='#38bdf8', marker='o', s=35, edgecolors='#0284c7', zorder=5, label=f'SPARC {gal_data["name"]}')
    ax1.plot(r_sim_kpc, rho_sim_norm, color='#f59e0b', marker='s', markersize=4.5, linewidth=2.2, zorder=6, label=f'Reotransductor Eon {engine.eon} Core')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(0.2, r_max_gal * 1.1)
    ax1.set_ylabel(r'Density $\rho(r)\ \ [M_\odot/\mathrm{pc}^3]$', color='#f8fafc', fontsize=11, fontweight='bold')
    ax1.set_title(f'A) Density Profile (Eon {engine.eon} vs SPARC NGC 2403)', color='#f8fafc', fontsize=12, fontweight='bold')
    leg1 = ax1.legend(loc='upper right', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=8.5)
    for t in leg1.get_texts(): t.set_color('#e2e8f0')

    # Panel 2: Logarithmic Density Slope
    ax2.axhline(0.0, color='#06b6d4', linestyle=':', linewidth=1.5, label='Flat Core ($\gamma = 0.0$)')
    ax2.axhline(-1.0, color='#ef4444', linestyle=':', linewidth=1.5, label='NFW Cusp ($\gamma = -1.0$)')
    ax2.axhline(-3.0, color='#94a3b8', linestyle=':', linewidth=1.2, label='Outer Envelope ($\gamma = -3.0$)')
    ax2.plot(r_theory, gamma_nfw_theory, color='#ef4444', linestyle='--', linewidth=1.8)
    ax2.plot(r_theory, gamma_core_theory, color='#06b6d4', linestyle='-.', linewidth=1.8)
    ax2.plot(r_sim_kpc, gamma_sim, color='#f59e0b', marker='o', markersize=4.5, linewidth=2.2, label=f'Simulation $\gamma_0 = {halo_metrics["inner_slope_gamma0"]:.3f}$')
    ax2.set_xscale('log')
    ax2.set_xlim(0.2, r_max_gal * 1.1)
    ax2.set_ylim(-3.5, 0.5)
    ax2.set_xlabel(r'Galactocentric Radius $r\ \ [\mathrm{kpc}]$', color='#f8fafc', fontsize=11, fontweight='bold')
    ax2.set_ylabel(r'Log Slope $\gamma(r) = \mathrm{d}\ln\rho/\mathrm{d}\ln r$', color='#f8fafc', fontsize=10, fontweight='bold')
    ax2.set_title('B) Inner Slope & Cusp-Core Resolution', color='#f8fafc', fontsize=12, fontweight='bold')
    leg2 = ax2.legend(loc='lower left', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=8.5)
    for t in leg2.get_texts(): t.set_color('#e2e8f0')

    # Panel 3: Universal Stacked Rotation Curve (175 SPARC Galaxies)
    r_norm = np.array(pop_stats["stacked_r_norm"])
    v_med = np.array(pop_stats["stacked_v_median"])
    v_p16 = np.array(pop_stats["stacked_v_p16"])
    v_p84 = np.array(pop_stats["stacked_v_p84"])

    ax3.fill_between(r_norm, v_p16, v_p84, color='#38bdf8', alpha=0.22, label=r'175 SPARC Galaxies ($1\sigma$ Dispersion)')
    ax3.plot(r_norm, v_med, color='#38bdf8', linestyle='--', linewidth=2.0, label='SPARC Catalog Median')
    # Simulation normalized rotation curve
    r_sim_norm = r_sim / max(1e-3, np.max(r_sim))
    v_sim_norm = v_sim / max(1e-3, np.max(v_sim))
    ax3.plot(r_sim_norm, v_sim_norm, color='#f59e0b', marker='s', markersize=4.5, linewidth=2.4, label='Reotransductor Universal Core')
    ax3.set_xlim(0.0, 1.05)
    ax3.set_ylim(0.0, 1.35)
    ax3.set_xlabel(r'Normalized Radius $r / R_{\mathrm{max}}$', color='#f8fafc', fontsize=11, fontweight='bold')
    ax3.set_ylabel(r'Normalized Velocity $V / V_{\mathrm{flat}}$', color='#f8fafc', fontsize=11, fontweight='bold')
    ax3.set_title(f'C) Universal Stacked Rotation Curves (N = {pop_stats["total_galaxies_evaluated"]})', color='#f8fafc', fontsize=12, fontweight='bold')
    leg3 = ax3.legend(loc='lower right', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=8.5)
    for t in leg3.get_texts(): t.set_color('#e2e8f0')

    # Panel 4: Population Goodness-of-Fit Delta Chi^2 Distribution
    delta_chi2 = np.array(pop_stats["delta_chi2_list"])
    bins = np.linspace(-15.0, 25.0, 30)
    ax4.hist(delta_chi2[delta_chi2 > 0], bins=bins, color='#06b6d4', alpha=0.85, edgecolor='#0891b2', label=f'Preferred Core: {pop_stats["core_preferred_count"]} ({pop_stats["core_preference_pct"]}%)')
    ax4.hist(delta_chi2[delta_chi2 <= 0], bins=bins, color='#ef4444', alpha=0.85, edgecolor='#dc2626', label=f'Preferred NFW: {pop_stats["nfw_preferred_count"]} ({pop_stats["nfw_preference_pct"]}%)')
    ax4.axvline(0.0, color='#f8fafc', linestyle='--', linewidth=1.5)
    ax4.set_xlabel(r'$\Delta\chi^2 = \chi^2_{\mathrm{NFW}} - \chi^2_{\mathrm{Core}}$ (Positive $\to$ Core Favored)', color='#f8fafc', fontsize=10, fontweight='bold')
    ax4.set_ylabel('Number of Galaxies', color='#f8fafc', fontsize=11, fontweight='bold')
    ax4.set_title('D) Model Selection Distribution (SPARC Catalog)', color='#f8fafc', fontsize=12, fontweight='bold')
    leg4 = ax4.legend(loc='upper right', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=8.5)
    for t in leg4.get_texts(): t.set_color('#e2e8f0')

    plt.suptitle(
        r'$\mathbf{SPARC\ 2020\ Population\ Benchmark\ (N=175)\ vs.\ Reotransductor\ Emergent\ Cored\ Halos}$',
        color='#f8fafc', fontsize=14, fontweight='bold', y=0.98
    )

    plt.savefig(output_fig, dpi=180, bbox_inches='tight', facecolor=fig.get_facecolor())
    pop_fig_path = "results/paper_II/halo_sparc_population_175.png"
    try:
        shutil.copyfile(output_fig, pop_fig_path)
    except Exception:
        pass
    plt.close(fig)

    print(f"• Output 4-Panel Population Figure: {output_fig} & {pop_fig_path}")
    print("=" * 75)
    return pop_stats


def process_all_existing_checkpoints_halo(checkpoints_dir: str = "checkpoints", galaxy_name: str = "NGC_2403"):
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
    os.makedirs(snapshots_dir, exist_ok=True)

    for npz_path in npz_files:
        try:
            engine = CosmologicalEngine.from_checkpoint(npz_path)

            print(f"\n--- Processing Halo: {os.path.basename(npz_path)} (Eon {engine.eon}, Grid {engine.grid_size}³) ---")
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
    parser.add_argument("--all-galaxies", action="store_true", help="Run full population benchmark across all 175 SPARC galaxies")
    parser.add_argument("--galaxy", type=str, default="NGC_2403", help="SPARC benchmark galaxy identifier (e.g. NGC_2403, DDO_154)")
    parser.add_argument("--output", type=str, default="results/paper_II/halo_cusp_core_sparc.png", help="Output PNG path")
    args = parser.parse_args()

    if args.all_galaxies:
        run_full_sparc_population_analysis(checkpoint_path=args.checkpoint_dir, output_fig=args.output)
    elif args.process_all:
        process_all_existing_checkpoints_halo(checkpoints_dir=args.checkpoint_dir, galaxy_name=args.galaxy)
    else:
        run_halo_comparison(
            checkpoint_path=args.checkpoint_dir,
            auto_resume=not args.from_scratch,
            steps=args.steps,
            galaxy_name=args.galaxy,
            output_fig=args.output
        )
