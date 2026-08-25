"""
Command-Line Utility and Publication Plot Generator:
Compares Reotransductor Cosmological Simulation with DESI 2024 & SDSS BOSS BAO Observations.
Generates results/paper_II/bao_comparison_desi.png and per-eon snapshot reports.
"""

import os
import glob
import shutil
import argparse
from typing import Dict, Any, Optional
import matplotlib.pyplot as plt
import numpy as np

from server.engine import CosmologicalEngine
from observational.bao_data import DESI2024BAOData
from observational.bao_analyzer import BAOSpatialCorrelationAnalyzer


def run_bao_comparison(
    checkpoint_path: str = "checkpoints",
    auto_resume: bool = True,
    steps: int = 0,
    output_fig: str = "results/paper_II/bao_comparison_desi.png"
) -> Dict[str, Any]:
    """Executes live BAO comparison on active cosmological state and outputs figure."""
    os.makedirs(os.path.dirname(output_fig) or ".", exist_ok=True)
    print("=" * 75)
    print("  🌌 REOTRANSDUCTOR OBSERVATIONAL PIPELINE: DESI 2024 BAO VALIDATION")
    print("=" * 75)
    target_checkpoint = None
    if os.path.isdir(checkpoint_path):
        bao_files = sorted(glob.glob(os.path.join(checkpoint_path, "bao_eon_*.npz")))
        if bao_files:
            target_checkpoint = bao_files[-1]
    elif os.path.isfile(checkpoint_path):
        target_checkpoint = checkpoint_path

    if target_checkpoint and os.path.exists(target_checkpoint):
        print(f"• Loading Target BAO/Cosmic Noon Checkpoint: '{target_checkpoint}'...")
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
    metrics = generate_eon_bao_report(engine, output_dir=snapshots_dir)
    print("\n" + "=" * 75)
    print("  📊 BAO SPATIAL CLUSTERING & STATISTICAL METRICS")
    print("=" * 75)
    print(f"• Correlation Length r_0:          {metrics['correlation_length_r0_mpc']:.2f} Mpc")
    print(f"• Detected Acoustic Peak r_BAO:    {metrics['bao_peak_radius_mpc']:.2f} Mpc")
    print(f"• DESI 2024 Goodness of Fit Chi^2: {metrics['chi2_stats']['chi2']} (Reduced Chi^2: {metrics['chi2_stats']['reduced_chi2']})")
    print(f"• Output Plot:                     {output_fig}")
    print("=" * 75)
    return metrics


def generate_eon_bao_report(engine, output_dir: str = "checkpoints/snapshots") -> Dict[str, Any]:
    """
    Generates complete BAO observational report for an active or checkpointed eon:
      1. Computes 3D spatial correlation function xi_rho(r) and memory xi_tau(r)
      2. Compares against official DESI 2024 DR1 and SDSS BOSS DR12 datasets
      3. Renders and saves publication plot bao_comparison_eon_N.png and results/paper_II/bao_comparison_desi.png
    """
    os.makedirs(output_dir, exist_ok=True)
    eon = engine.eon
    scale_factor = float(engine.scale_factor)

    # 1. Ingest Observational Datasets
    desi_data = DESI2024BAOData()
    r_desi, xi_desi, err_desi = desi_data.get_desi_dataset()
    r_boss, xi_boss = desi_data.get_boss_dataset()
    r_lcdm, xi_lcdm = desi_data.get_theoretical_lcdm_correlation()

    # 2. Compute 3D Spatial Correlation from Simulation Grid
    analyzer = BAOSpatialCorrelationAnalyzer(
        grid_size=engine.grid_size,
        box_size_mpc=getattr(engine, 'box_size_mpc', 100.0),
        n_bins=32
    )
    sim_metrics = analyzer.evaluate_cosmological_fields(
        rho_3d=engine.rho,
        tau_3d=engine.tau,
        scale_factor=scale_factor
    )

    r_sim = np.array(sim_metrics["r_comoving_mpc"])
    xi_sim = np.array(sim_metrics["xi_rho"])

    # True comoving separation distance within the simulation box (0 to L_box / 2)
    r_sim_plot = r_sim
    xi_sim_plot = xi_sim

    # 3. Compute Chi^2 Goodness of Fit against DESI dataset overlapping bins
    chi2_stats = desi_data.compute_chi2(r_sim_plot, xi_sim_plot)

    # 4. Render Publication Figure
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 8),
        gridspec_kw={'height_ratios': [2.2, 1.0], 'hspace': 0.28}
    )
    fig.patch.set_facecolor('#0d131f')

    for ax in (ax1, ax2):
        ax.set_facecolor('#131b2e')
        ax.grid(True, linestyle=':', alpha=0.35, color='#475569')
        ax.tick_params(colors='#94a3b8', labelsize=10)
        for spine in ax.spines.values():
            spine.set_edgecolor('#334155')

    # -------------------------------------------------------------------------
    # Panel 1: Monopole Two-Point Correlation r^2 * xi(r)
    # -------------------------------------------------------------------------
    # LCDM Linear Theory Baseline
    ax1.plot(
        r_lcdm, (r_lcdm**2) * xi_lcdm,
        color='#94a3b8', linestyle='--', linewidth=1.8,
        label=r'Standard $\Lambda$CDM Baseline ($r_s = 103.0\ h^{-1}\mathrm{Mpc}$)'
    )

    # SDSS BOSS DR12 Data
    ax1.scatter(
        r_boss, (r_boss**2) * xi_boss,
        color='#3b82f6', marker='D', s=35, alpha=0.8,
        label='SDSS BOSS DR12 Consensus (Alam et al.)'
    )

    # DESI 2024 DR1 Data with 1-sigma Error Bars
    y_desi = (r_desi**2) * xi_desi
    y_err_desi = (r_desi**2) * err_desi
    ax1.errorbar(
        r_desi, y_desi, yerr=y_err_desi,
        fmt='o', color='#06b6d4', ecolor='#22d3ee', elinewidth=1.6, capsize=3.5,
        markersize=6, label='DESI 2024 DR1 LRG/QSO (Official Points)'
    )

    # Reotransductor Simulation Curve
    r2_xi_sim = (r_sim_plot**2) * xi_sim_plot
    ax1.plot(
        r_sim_plot, r2_xi_sim,
        color='#f59e0b', marker='s', markersize=4.5, linewidth=2.4,
        label=f'Reotransductor Eon {eon} Simulation ($L_{{\\mathrm{{box}}}} = 100\\ h^{{-1}}\\mathrm{{Mpc}}$)'
    )

    ax1.set_xlim(0.0, 155.0)
    ax1.set_ylabel(r'$r^2 \, \xi(r)\ \ [h^{-2}\mathrm{Mpc}^2]$', color='#f8fafc', fontsize=12, fontweight='bold')
    ax1.set_title(
        f'Baryon Acoustic Oscillations (BAO) Spatial Correlation — Eon {eon} vs. DESI 2024 DR1\n'
        f'(Scale Factor a = {scale_factor:.3f} | Exploratory Lattice Benchmark)',
        color='#f8fafc', fontsize=12, fontweight='bold', pad=12
    )
    legend = ax1.legend(loc='upper right', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=9.5)
    for text in legend.get_texts():
        text.set_color('#e2e8f0')

    # -------------------------------------------------------------------------
    # Panel 2: Simulation Matter Correlation xi_sim(r)
    # -------------------------------------------------------------------------
    ax2.plot(
        r_sim_plot, xi_sim_plot,
        color='#f59e0b', marker='o', markersize=4, linewidth=1.8,
        label=r'Simulation $\xi_{\rho}(r)$ (Unscaled Comoving Frame)'
    )
    ax2.axhline(0.0, color='#94a3b8', linestyle='--', linewidth=1.2)

    ax2.set_xlim(0.0, 155.0)
    ax2.set_xlabel(r'Comoving Separation $r\ \ [h^{-1}\mathrm{Mpc}]$', color='#f8fafc', fontsize=11, fontweight='bold')
    ax2.set_ylabel(r'$\xi(r)$', color='#f8fafc', fontsize=10, fontweight='bold')

    leg2 = ax2.legend(loc='upper right', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=8.5)
    for text in leg2.get_texts():
        text.set_color('#e2e8f0')

    # Save Eon Plot
    eon_bao_path = os.path.join(output_dir, f"bao_comparison_eon_{eon}.png")
    plt.savefig(eon_bao_path, dpi=180, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)

    # Save Latest Preview in results/paper_II/
    latest_bao_path = "results/paper_II/bao_comparison_desi.png"
    try:
        os.makedirs("assets", exist_ok=True)
        shutil.copyfile(eon_bao_path, latest_bao_path)
    except Exception:
        pass

    return {
        "eon": eon,
        "scale_factor": scale_factor,
        "correlation_length_r0_mpc": sim_metrics["correlation_length_r0_mpc"],
        "bao_peak_radius_mpc": sim_metrics["bao_peak_radius_mpc"],
        "bao_peak_amplitude": sim_metrics["bao_peak_amplitude"],
        "chi2_stats": chi2_stats,
        "bao_comparison_path": eon_bao_path
    }


def process_all_existing_checkpoints_bao(checkpoints_dir: str = "checkpoints"):
    """Processes all historical BAO checkpoints (bao_eon_*.npz or eon_*.npz) and generates complete plots."""
    # Prioritize dedicated BAO epoch checkpoints (a ~ 2.0)
    npz_files = sorted(glob.glob(os.path.join(checkpoints_dir, "bao_eon_*.npz")))
    if not npz_files:
        npz_files = sorted(glob.glob(os.path.join(checkpoints_dir, "eon_*.npz")))

    if not npz_files:
        print(f"No BAO or eon checkpoints found in '{checkpoints_dir}'.")
        return

    print(f"• Found {len(npz_files)} BAO epoch checkpoints in '{checkpoints_dir}'. Processing BAO reports...")
    snapshots_dir = os.path.join(checkpoints_dir, "snapshots")

    for npz_path in npz_files:
        try:
            engine = CosmologicalEngine.from_checkpoint(npz_path)

            print(f"\n--- Processing BAO: {os.path.basename(npz_path)} (Eon {engine.eon}) ---")
            metrics = generate_eon_bao_report(engine, output_dir=snapshots_dir)
            print(f"  * Peak r_BAO: {metrics['bao_peak_radius_mpc']} Mpc | Chi2: {metrics['chi2_stats']['chi2']} (Reduced: {metrics['chi2_stats']['reduced_chi2']})")
        except Exception as ex:
            print(f"  * Error processing {npz_path}: {ex}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Reotransductor vs. DESI 2024 BAO Observational Validation")
    parser.add_argument("--checkpoint-dir", type=str, default="checkpoints", help="Directory with latest checkpoint or snapshot")
    parser.add_argument("--steps", type=int, default=0, help="Additional integration steps to simulate before evaluation")
    parser.add_argument("--from-scratch", action="store_true", help="Start from blank initial conditions rather than active checkpoint")
    parser.add_argument("--process-all", action="store_true", help="Process and generate BAO plots for all historical eon_*.npz checkpoints")
    parser.add_argument("--output", type=str, default="results/paper_II/bao_comparison_desi.png", help="Output PNG path")
    args = parser.parse_args()

    if args.process_all:
        process_all_existing_checkpoints_bao(checkpoints_dir=args.checkpoint_dir)
    else:
        run_bao_comparison(
            checkpoint_path=args.checkpoint_dir,
            auto_resume=not args.from_scratch,
            steps=args.steps,
            output_fig=args.output
        )
