"""
Command-Line Utility and Publication Plot Generator:
Compares Reotransductor Cosmological Simulation with Pantheon+ (2022) Supernovae Catalog
and evaluates the resolution of the Cosmological Hubble Tension (assets/pantheon_hubble_tension.png).
"""

import os
import glob
import shutil
import argparse
from typing import Dict, Any, Optional
import matplotlib.pyplot as plt
import numpy as np

from server.engine import CosmologicalEngine
from observational.pantheon_data import PantheonSupernovaeData
from observational.hubble_tension import HubbleTensionAnalyzer


def run_pantheon_comparison(
    checkpoint_path: str = "checkpoints",
    auto_resume: bool = True,
    steps: int = 0,
    mode: str = "full",
    output_fig: str = "assets/pantheon_hubble_tension.png"
) -> Dict[str, Any]:
    """Executes live Hubble tension and Pantheon+ supernovae comparison on active cosmological state."""
    os.makedirs(os.path.dirname(output_fig) or ".", exist_ok=True)
    print("=" * 75)
    print(f"  🌌 REOTRANSDUCTOR OBSERVATIONAL PIPELINE: HUBBLE TENSION & PANTHEON+ SNe ({mode.upper()})")
    print("=" * 75)
    # Prioritize dedicated local universe epoch checkpoints (a ~ 4.5)
    target_checkpoint = None
    if os.path.isdir(checkpoint_path):
        pantheon_files = sorted(glob.glob(os.path.join(checkpoint_path, "pantheon_eon_*.npz")))
        if pantheon_files:
            target_checkpoint = pantheon_files[-1]
    elif os.path.isfile(checkpoint_path):
        target_checkpoint = checkpoint_path

    if target_checkpoint and os.path.exists(target_checkpoint):
        print(f"• Loading Target Local Universe Checkpoint: '{target_checkpoint}'...")
        data = np.load(target_checkpoint)
        grid_n = int(data['rho'].shape[0])
        engine = CosmologicalEngine(grid_size=grid_n, checkpoint_dir=os.path.dirname(target_checkpoint), auto_resume=False)
        engine.rho = data['rho']
        engine.v_x = data.get('v_x', None)
        engine.v_y = data.get('v_y', None)
        engine.v_z = data.get('v_z', None)
        engine.T = data.get('T', None)
        engine.I = data.get('I', None)
        engine.tau = data['tau']
        if 'phi' in data:
            engine.phi = data['phi']
        engine.scale_factor = float(data['scale_factor'])
        engine.eon = int(data['eon'])
        engine.total_steps = int(data.get('total_steps', 0))
    else:
        print(f"• Loading Cosmological State (checkpoint_dir='{checkpoint_path}', auto_resume={auto_resume})...")
        engine = CosmologicalEngine(checkpoint_dir=checkpoint_path, auto_resume=auto_resume)
        if steps > 0:
            print(f"• Evolving additional {steps} integration steps...")
            for _ in range(steps):
                engine.step()

    print(f"• Evaluated State: Eon {engine.eon} | Steps: {engine.total_steps:,} | Scale Factor a = {engine.scale_factor:.3f}")
    snapshots_dir = "checkpoints/snapshots"
    metrics = generate_eon_pantheon_report(engine, mode=mode, output_dir=snapshots_dir)
    print("\n" + "=" * 75)
    print("  📊 HUBBLE TENSION & ENVIRONMENTAL EXPANSION DIAGNOSTICS")
    print("=" * 75)
    print(f"• Pantheon+ Sample Size:            {metrics['total_sne']} Supernovae Type Ia")
    print(f"• Planck 2018 Baseline (CMB/Void): {metrics['h0_background_planck']} km/s/Mpc")
    print(f"• SH0ES 2022 Observed (Cluster):   {metrics['h0_observed_shoes']} km/s/Mpc (Tension: +{metrics['delta_h0_observed']} km/s/Mpc)")
    print(f"• Reotransductor Predicted H_0:     {metrics['h0_predicted_reotransductor']} km/s/Mpc (Delta: +{metrics['delta_h0_predicted']} km/s/Mpc)")
    print(f"• Tension Resolution Percentage:    {metrics['tension_resolution_pct']}% ({'✅ RESOLVED' if metrics['is_tension_mitigated'] else '⚠️ PARTIAL'})")
    print(f"• Environmental Gradient (dH0/dlog): {metrics['environmental_field']['environmental_gradient_slope']} km/s/Mpc/dex")
    print(f"• Chi^2 Goodness of Fit vs SNe:     {metrics['chi2_reotransductor']} (Reotransductor) vs {metrics['chi2_planck_static']} (Planck Static)")
    print(f"• Output Plot:                      {output_fig}")
    print("=" * 75)
    return metrics


def generate_eon_pantheon_report(
    engine,
    mode: str = "full",
    output_dir: str = "checkpoints/snapshots"
) -> Dict[str, Any]:
    """
    Generates complete Hubble Tension and Pantheon+ Supernovae observational report:
      1. Evaluates local dilated expansion rate H_0(x) across cosmic web environments
      2. Ingests official Pantheon+ (2022) distance modulus data (1,701 full sample or binned)
      3. Computes Chi^2 and residual curves
      4. Renders and saves 3-panel publication figure
    """
    os.makedirs(output_dir, exist_ok=True)
    eon = engine.eon
    scale_factor = float(engine.scale_factor)

    # 1. Ingest Pantheon+ Observational Data
    pantheon = PantheonSupernovaeData(mode=mode)
    data = pantheon.get_dataset()
    z_sne = data["z"]
    mu_obs = data["mu_obs"]
    err_mu = data["err_mu"]

    # 2. Evaluate Engine Hubble Rate & Environmental Field
    analyzer = HubbleTensionAnalyzer(
        h0_bg=pantheon.benchmarks["h0_planck"],
        h0_local_obs=pantheon.benchmarks["h0_shoes"]
    )
    hubble_metrics = analyzer.evaluate_engine_state(engine)
    env_field = hubble_metrics["environmental_field"]

    h0_pred = hubble_metrics["h0_predicted_reotransductor"]
    h0_planck = hubble_metrics["h0_background_planck"]
    h0_shoes = hubble_metrics["h0_observed_shoes"]

    # 3. Compute Theoretical Distance Modulus Curves
    z_theory = np.linspace(0.005, 1.4, 150)
    mu_planck_theory = pantheon.distance_modulus(z_theory, h0=h0_planck)
    mu_shoes_theory = pantheon.distance_modulus(z_theory, h0=h0_shoes)

    # Reotransductor environmental transition curve:
    # At low z (local galaxies), SNe are situated inside dense host halos (H0 -> H0_local)
    # At high z, light traverses vast cosmic voids and background Hubble flow (H0 -> H0_CMB)
    w_local = np.exp(-z_theory / 0.15)  # Local clustering transition kernel
    h0_effective_z = h0_planck + (h0_pred - h0_planck) * w_local
    
    mu_reotransductor_theory = np.zeros_like(z_theory)
    for i, z in enumerate(z_theory):
        mu_reotransductor_theory[i] = pantheon.distance_modulus(np.array([z]), h0=h0_effective_z[i])[0]

    # Compute model predictions at discrete supernovae redshifts
    w_sne = np.exp(-z_sne / 0.15)
    h0_sne_eff = h0_planck + (h0_pred - h0_planck) * w_sne
    mu_reotransductor_sne = np.zeros_like(z_sne)
    mu_planck_sne = pantheon.distance_modulus(z_sne, h0=h0_planck)
    mu_shoes_sne = pantheon.distance_modulus(z_sne, h0=h0_shoes)

    for i, z in enumerate(z_sne):
        mu_reotransductor_sne[i] = pantheon.distance_modulus(np.array([z]), h0=h0_sne_eff[i])[0]

    # Chi^2 Goodness of Fit
    chi2_reotransductor = round(pantheon.compute_chi2(mu_reotransductor_sne, mu_obs, err_mu), 2)
    chi2_planck = round(pantheon.compute_chi2(mu_planck_sne, mu_obs, err_mu), 2)
    chi2_shoes = round(pantheon.compute_chi2(mu_shoes_sne, mu_obs, err_mu), 2)

    residuals_reotransductor = (mu_obs - mu_reotransductor_sne) / err_mu

    # 4. Render 3-Panel Publication Figure
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(10, 11),
        gridspec_kw={'height_ratios': [1.8, 1.3, 1.2], 'hspace': 0.28}
    )
    fig.patch.set_facecolor('#0d131f')

    for ax in (ax1, ax2, ax3):
        ax.set_facecolor('#131b2e')
        ax.grid(True, linestyle=':', alpha=0.35, color='#475569')
        ax.tick_params(colors='#94a3b8', labelsize=10)
        for spine in ax.spines.values():
            spine.set_edgecolor('#334155')

    # -------------------------------------------------------------------------
    # Panel 1: Hubble Diagram mu(z)
    # -------------------------------------------------------------------------
    ax1.plot(
        z_theory, mu_planck_theory,
        color='#ef4444', linestyle='--', linewidth=1.8,
        label=f'Planck 2018 Baseline ($H_0 = {h0_planck:.2f}\\ \\mathrm{{km/s/Mpc}}$, $\\chi^2 = {chi2_planck}$)'
    )
    ax1.plot(
        z_theory, mu_shoes_theory,
        color='#06b6d4', linestyle='-.', linewidth=1.8,
        label=f'SH0ES 2022 Local Ladder ($H_0 = {h0_shoes:.2f}\\ \\mathrm{{km/s/Mpc}}$, $\\chi^2 = {chi2_shoes}$)'
    )
    ax1.plot(
        z_theory, mu_reotransductor_theory,
        color='#f59e0b', linewidth=2.4, zorder=5,
        label=f'Reotransductor Environmental Model ($H_0(\\mathbf{{x}})$, $\\chi^2 = {chi2_reotransductor}$)'
    )
    if len(z_sne) > 100:
        ax1.errorbar(
            z_sne, mu_obs, yerr=err_mu,
            fmt='o', color='#38bdf8', ecolor='#0284c7', elinewidth=0.8, alpha=0.35, capsize=0,
            markersize=3.0, zorder=4, label=f'Pantheon+ (2022) Full Sample ($N = {len(z_sne):,}$ SNe Ia)'
        )
    else:
        ax1.errorbar(
            z_sne, mu_obs, yerr=err_mu,
            fmt='o', color='#38bdf8', ecolor='#0284c7', elinewidth=1.6, capsize=3.5,
            markersize=5.5, zorder=6, label=f'Pantheon+ (2022) Calibration Subset ($N = {len(z_sne)}$ Bins)'
        )

    ax1.set_xscale('log')
    ax1.set_xlim(0.008, 1.5)
    ax1.set_ylim(32.5, 46.0)
    ax1.set_ylabel(r'Distance Modulus $\mu(z)$ [mag]', color='#f8fafc', fontsize=11, fontweight='bold')
    ax1.set_title(
        f'Environmental Hubble Rate Analysis — Eon {eon} vs. Pantheon+ SNe Ia (2022)\n'
        f'(Cluster $H_0 = {h0_pred:.2f}\\ \\mathrm{{km/s/Mpc}}$ | Sample: {len(z_sne):,} SNe Ia)',
        color='#f8fafc', fontsize=12, fontweight='bold', pad=10
    )
    leg1 = ax1.legend(loc='lower right', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=9)
    for t in leg1.get_texts():
        t.set_color('#e2e8f0')

    # -------------------------------------------------------------------------
    # Panel 2: Environmental H_0 Gradient vs Local Density Contrast
    # -------------------------------------------------------------------------
    ax2.axhline(h0_planck, color='#ef4444', linestyle='--', linewidth=1.5, label=f'Planck CMB Void Baseline ({h0_planck} km/s/Mpc)')
    ax2.axhline(h0_shoes, color='#06b6d4', linestyle='-.', linewidth=1.5, label=f'SH0ES Local Cluster Target ({h0_shoes} km/s/Mpc)')

    centers = np.array(env_field["density_bin_centers"])
    h0_binned = np.array(env_field["binned_h0"])
    h0_errs = np.array(env_field["binned_h0_err"])

    ax2.errorbar(
        centers, h0_binned, yerr=h0_errs,
        fmt='s', color='#f59e0b', ecolor='#d97706', elinewidth=1.8, capsize=3.5,
        markersize=6, label=f'Simulation Environmental $H_0(\\delta)$ (Slope = {env_field["environmental_gradient_slope"]:.2f})'
    )

    # Shade environment domains
    ax2.axvspan(-1.0, -0.3, alpha=0.10, color='#3b82f6', label='Cosmic Voids')
    ax2.axvspan(-0.3, 0.5, alpha=0.10, color='#10b981', label='Filaments & Walls')
    ax2.axvspan(0.5, 1.5, alpha=0.10, color='#f59e0b', label='Clusters & Halos')

    ax2.set_xlim(-1.0, 1.5)
    ax2.set_ylim(64.0, 76.0)
    ax2.set_xlabel(r'Local Density Contrast $\log_{10}(\rho / \bar{\rho})$', color='#f8fafc', fontsize=11, fontweight='bold')
    ax2.set_ylabel(r'$H_0(\mathbf{x})\ \ [\mathrm{km/s/Mpc}]$', color='#f8fafc', fontsize=10, fontweight='bold')
    leg2 = ax2.legend(loc='lower right', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=8.5, ncol=2)
    for t in leg2.get_texts():
        t.set_color('#e2e8f0')

    # -------------------------------------------------------------------------
    # Panel 3: Standardized Residuals Delta mu / sigma_mu
    # -------------------------------------------------------------------------
    ax3.axhline(0.0, color='#94a3b8', linestyle='-', linewidth=1.2)
    ax3.axhspan(-1.0, 1.0, alpha=0.15, color='#10b981', label=r'$\pm 1\sigma$ Observational Confidence')
    ax3.axhspan(-2.0, 2.0, alpha=0.08, color='#3b82f6', label=r'$\pm 2\sigma$ Confidence')

    pt_size = 12 if len(z_sne) > 100 else 45
    pt_alpha = 0.40 if len(z_sne) > 100 else 0.85
    ax3.scatter(
        z_sne, residuals_reotransductor,
        color='#38bdf8', marker='o', s=pt_size, alpha=pt_alpha, edgecolors='none' if len(z_sne) > 100 else '#0284c7', zorder=5,
        label=r'Normalized Residuals $(\mu_{\mathrm{obs}} - \mu_{\mathrm{model}}) / \sigma_\mu$'
    )

    ax3.set_xscale('log')
    ax3.set_xlim(0.008, 1.5)
    ax3.set_ylim(-3.0, 3.0)
    ax3.set_xlabel(r'Redshift $z_{\mathrm{CMB}}$', color='#f8fafc', fontsize=11, fontweight='bold')
    ax3.set_ylabel(r'Residual $[\sigma]$', color='#f8fafc', fontsize=10, fontweight='bold')
    leg3 = ax3.legend(loc='lower left', framealpha=0.85, facecolor='#0b1120', edgecolor='#334155', fontsize=8.5)
    for t in leg3.get_texts():
        t.set_color('#e2e8f0')

    # Save Eon Plot
    eon_pantheon_path = os.path.join(output_dir, f"pantheon_comparison_eon_{eon}.png")
    plt.savefig(eon_pantheon_path, dpi=180, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close(fig)

    # Copy to assets/
    latest_pantheon_path = "assets/pantheon_hubble_tension.png"
    try:
        os.makedirs("assets", exist_ok=True)
        shutil.copyfile(eon_pantheon_path, latest_pantheon_path)
    except Exception:
        pass

    return {
        "eon": eon,
        "scale_factor": scale_factor,
        "total_sne": len(z_sne),
        "h0_background_planck": h0_planck,
        "h0_observed_shoes": h0_shoes,
        "h0_predicted_reotransductor": h0_pred,
        "delta_h0_observed": hubble_metrics["delta_h0_observed"],
        "delta_h0_predicted": hubble_metrics["delta_h0_predicted"],
        "tension_resolution_pct": hubble_metrics["tension_resolution_pct"],
        "is_tension_mitigated": hubble_metrics["is_tension_mitigated"],
        "environmental_field": env_field,
        "chi2_reotransductor": chi2_reotransductor,
        "chi2_planck_static": chi2_planck,
        "chi2_shoes_static": chi2_shoes,
        "pantheon_comparison_path": eon_pantheon_path
    }


def process_all_existing_checkpoints_pantheon(checkpoints_dir: str = "checkpoints", mode: str = "full"):
    """Processes all historical Pantheon/Universe Local checkpoints (pantheon_eon_*.npz or eon_*.npz)."""
    # Prioritize dedicated local universe epoch checkpoints (a ~ 4.5)
    npz_files = sorted(glob.glob(os.path.join(checkpoints_dir, "pantheon_eon_*.npz")))
    if not npz_files:
        npz_files = sorted(glob.glob(os.path.join(checkpoints_dir, "eon_*.npz")))

    if not npz_files:
        print(f"No Pantheon or eon checkpoints found in '{checkpoints_dir}'.")
        return

    print(f"• Found {len(npz_files)} Pantheon checkpoints in '{checkpoints_dir}'. Processing Hubble Tension reports...")
    snapshots_dir = os.path.join(checkpoints_dir, "snapshots")
    engine = CosmologicalEngine(checkpoint_dir=checkpoints_dir, auto_resume=False)

    for npz_path in npz_files:
        try:
            data = np.load(npz_path)
            engine.grid_size = int(data['rho'].shape[0])
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

            print(f"\n--- Processing Hubble: {os.path.basename(npz_path)} (Eon {engine.eon}, Grid {engine.grid_size}³) ---")
            metrics = generate_eon_pantheon_report(engine, mode=mode, output_dir=snapshots_dir)
            print(f"  * H0 Predicted: {metrics['h0_predicted_reotransductor']} km/s/Mpc | Resolution: {metrics['tension_resolution_pct']}% | Chi2: {metrics['chi2_reotransductor']}")
        except Exception as ex:
            print(f"  * Error processing {npz_path}: {ex}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Reotransductor vs. Pantheon+ (2022) Hubble Tension Validation")
    parser.add_argument("--checkpoint-dir", type=str, default="checkpoints", help="Directory with checkpoints")
    parser.add_argument("--steps", type=int, default=0, help="Additional steps to simulate before evaluation")
    parser.add_argument("--from-scratch", action="store_true", help="Start from blank initial conditions")
    parser.add_argument("--process-all", action="store_true", help="Process and generate Pantheon plots for all historical checkpoints")
    parser.add_argument("--binned", action="store_true", help="Use 13-bin calibration subset instead of full 1,701 SNe sample")
    parser.add_argument("--output", type=str, default="assets/pantheon_hubble_tension.png", help="Output PNG path")
    args = parser.parse_args()

    mode = "binned" if args.binned else "full"
    if args.process_all:
        process_all_existing_checkpoints_pantheon(checkpoints_dir=args.checkpoint_dir, mode=mode)
    else:
        run_pantheon_comparison(
            checkpoint_path=args.checkpoint_dir,
            auto_resume=not args.from_scratch,
            steps=args.steps,
            mode=mode,
            output_fig=args.output
        )
