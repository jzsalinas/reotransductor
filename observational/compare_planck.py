"""
Command-Line Utility and Publication Plot Generator:
Compares Reotransductor Cosmological Simulation with ESA Planck 2018 CMB Data.
Generates assets/planck_comparison_spectrum.png.
"""

import os
import shutil
import argparse
from typing import Dict, Any
import matplotlib.pyplot as plt
import numpy as np
from server.engine import CosmologicalEngine
from observational.planck_data import Planck2018Data
from observational.cmb_analyzer import CMBSphericalHarmonicsAnalyzer
from observational.hubble_tension import HubbleTensionAnalyzer


def run_planck_comparison(
    checkpoint_path: str = "checkpoints",
    auto_resume: bool = True,
    steps: int = 0,
    output_fig: str = "assets/planck_comparison_spectrum.png"
):
    """Executes live comparison on active state and outputs figure."""
    os.makedirs(os.path.dirname(output_fig), exist_ok=True)
    print("=" * 70)
    print("  🌌 REOTRANSDUCTOR OBSERVATIONAL PIPELINE: ESA PLANCK 2018 VALIDATION")
    print("=" * 70)
    print(f"• Loading Cosmological State (checkpoint_dir='{checkpoint_path}', auto_resume={auto_resume})...")
    engine = CosmologicalEngine(checkpoint_dir=checkpoint_path, auto_resume=auto_resume)
    if steps > 0:
        print(f"• Evolving additional {steps} integration steps...")
        for _ in range(steps):
            engine.step()

    print(f"• Evaluated State: Eon {engine.eon} | Total Steps: {engine.total_steps:,} | Scale Factor a = {engine.scale_factor:.3f} | Fossil Proper Time tau_max = {np.max(engine.tau):.1f}")
    metrics = generate_eon_observational_report(engine, output_dir=os.path.dirname(output_fig) or ".")
    print("\n" + "=" * 70)
    print("  📊 OBSERVATIONAL COMPARISON & STATISTICAL METRICS")
    print("=" * 70)
    print(f"• Quadrupole (ell=2) Power:        {metrics['quadrupole_C2']:.4f}")
    print(f"• Octopole (ell=3) Power:          {metrics['octopole_C3']:.4f}")
    print(f"• Quadrupole/Octopole Ratio:       {metrics['ratio_C2_C3']:.3f} (Low-ell suppression: {metrics['is_quadrupole_suppressed']})")
    print(f"• Planck Goodness of Fit Chi^2:    {metrics['planck_chi2']} (Reduced Chi^2: {metrics['planck_reduced_chi2']})")
    print(f"• Predicted Local H_0:             {metrics['h0_predicted']:.2f} km/s/Mpc")
    print(f"• Hubble Tension Resolution:       {metrics['h0_tension_resolution_pct']:.1f}%")
    print("=" * 70)
    return metrics


def generate_eon_observational_report(engine, output_dir: str = "checkpoints/snapshots") -> Dict[str, Any]:
    """
    Generates complete observational report for an active or checkpointed eon:
      1. CMB Mollweide projection map (cmb_mollweide_eon_N.png)
      2. ESA Planck 2018 comparative power spectrum & Hubble tension (planck_comparison_eon_N.png)
      3. Statistical metrics dictionary for history logging.
    """
    os.makedirs(output_dir, exist_ok=True)
    eon = engine.eon
    time_myr = float(engine.units.time_code_to_myr(np.max(engine.tau))) if hasattr(engine, 'units') else 0.0

    # 1. Ingest Planck 2018 Reference Data
    planck = Planck2018Data()
    ell_obs, D_obs, err_obs = planck.get_binned_spectrum()
    ell_lcdm, D_lcdm = planck.get_theoretical_lcdm_spectrum(ell_max=100)

    # 2. Extract Celestial Sphere and Spherical Harmonics Decomposition
    analyzer = CMBSphericalHarmonicsAnalyzer(n_theta=48, n_phi=96, ell_max=24)
    celestial_map = analyzer.extract_celestial_sphere_from_grid(
        tau_3d=engine.tau,
        T_3d=engine.T,
        rho_3d=engine.rho,
        v_x=engine.v_x,
        v_y=engine.v_y,
        v_z=engine.v_z,
        c_light=engine.C_LIGHT
    )

    spectrum = analyzer.compute_angular_power_spectrum(celestial_map, temperature_scale_uK=27.255)
    ell_sim = spectrum["ell"]
    D_sim = spectrum["D_ell"]

    # 3. Evaluate Hubble Tension
    ht_analyzer = HubbleTensionAnalyzer()
    ht_metrics = ht_analyzer.evaluate_engine_state(engine)
    chi2_stats = planck.compute_chi2(ell_sim, D_sim)

    # 4. Render and Save CMB Mollweide Map
    mollweide_eon_path = os.path.join(output_dir, f"cmb_mollweide_eon_{eon}.png")
    analyzer.render_mollweide_plot(
        celestial_map=celestial_map,
        output_path=mollweide_eon_path,
        eon=eon,
        scale_factor=float(engine.scale_factor),
        time_myr=time_myr
    )
    # Also save as latest preview in assets/
    latest_mollweide_path = "assets/latest_cmb_mollweide.png"
    try:
        os.makedirs("assets", exist_ok=True)
        shutil.copyfile(mollweide_eon_path, latest_mollweide_path)
    except Exception:
        pass

    # 5. Render and Save Planck Power Spectrum & Hubble Tension Plot
    planck_comp_path = os.path.join(output_dir, f"planck_comparison_eon_{eon}.png")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=False, gridspec_kw={'height_ratios': [2, 1.2], 'hspace': 0.32})
    fig.patch.set_facecolor('#0f172a')
    ax1.set_facecolor('#1e293b')
    ax2.set_facecolor('#1e293b')

    # Panel 1: Angular Power Spectrum D_ell
    ax1.errorbar(ell_obs[ell_obs <= 80], D_obs[ell_obs <= 80], yerr=err_obs[ell_obs <= 80], fmt='o', color='#38bdf8',
                 ecolor='#0284c7', elinewidth=1.5, capsize=3, label='ESA Planck 2018 Legacy Data (PR3)', zorder=4)
    ax1.plot(ell_lcdm[ell_lcdm <= 80], D_lcdm[ell_lcdm <= 80], color='#94a3b8', linestyle='--', linewidth=1.5,
             label=r'Standard $\Lambda$CDM Baseline ($H_0=67.36$)', zorder=2)
    ax1.plot(ell_sim, D_sim, color='#f59e0b', marker='s', markersize=5, linewidth=2,
             label=f'Reotransductor Eon {eon} Prediction ($D_\\ell$)', zorder=5)

    ax1.set_title(f'CMB Angular Power Spectrum Comparison — Eon {eon} vs. ESA Planck 2018', color='#f8fafc', fontsize=12, fontweight='bold', pad=10)
    ax1.set_ylabel(r'$D_\ell = \frac{\ell(\ell+1)}{2\pi} C_\ell \; [\mu\mathrm{K}^2]$', color='#f8fafc', fontsize=10)
    ax1.tick_params(colors='#cbd5e1', labelsize=9)
    ax1.grid(True, linestyle=':', alpha=0.3, color='#64748b')
    ax1.legend(loc='upper left', facecolor='#0f172a', edgecolor='#334155', labelcolor='#f8fafc', fontsize=9)

    # Panel 2: Hubble Tension
    categories = ['Planck 2018 (CMB)', f'Reotransductor (Eon {eon})', 'SH0ES 2022 (Local)']
    h0_values = [ht_metrics['h0_background_planck'], ht_metrics['h0_predicted_reotransductor'], ht_metrics['h0_observed_shoes']]
    h0_errors = [0.54, 0.45, 1.04]
    colors = ['#38bdf8', '#f59e0b', '#10b981']

    bars = ax2.barh(categories, h0_values, xerr=h0_errors, color=colors, alpha=0.85, capsize=4, edgecolor='#ffffff', height=0.5)
    ax2.set_xlim(60.0, 85.0)
    ax2.set_xlabel(r'Hubble Expansion Rate $H_0 \; [\mathrm{km/s/Mpc}]$', color='#f8fafc', fontsize=10)
    ax2.set_title(r'Hubble Tension Resolution via Dissipative Emergent Time $\Delta\tau = \kappa_0 \int \sigma \, dt$', color='#f8fafc', fontsize=11, pad=8)
    ax2.tick_params(colors='#cbd5e1', labelsize=9)
    ax2.grid(True, linestyle=':', alpha=0.3, color='#64748b', axis='x')

    for bar in bars:
        w = bar.get_width()
        ax2.text(w + 1.0, bar.get_y() + bar.get_height() / 2, f'{w:.2f} km/s/Mpc', color='#f8fafc', va='center', fontsize=9, fontweight='bold')

    plt.savefig(planck_comp_path, dpi=180, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close()
    print(f"• Planck comparison plot successfully saved to: {planck_comp_path}")

    # Also copy to assets/planck_comparison_spectrum.png
    latest_spectrum_path = "assets/planck_comparison_spectrum.png"
    try:
        shutil.copyfile(planck_comp_path, latest_spectrum_path)
    except Exception:
        pass

    return {
        "quadrupole_C2": round(float(spectrum["quadrupole_C2"]), 4),
        "octopole_C3": round(float(spectrum["octopole_C3"]), 4),
        "ratio_C2_C3": spectrum["ratio_C2_C3"],
        "is_quadrupole_suppressed": spectrum["is_quadrupole_suppressed"],
        "planck_chi2": chi2_stats["chi2"],
        "planck_reduced_chi2": chi2_stats["reduced_chi2"],
        "h0_predicted": ht_metrics["h0_predicted_reotransductor"],
        "h0_tension_resolution_pct": ht_metrics["tension_resolution_pct"],
        "cmb_mollweide_path": mollweide_eon_path,
        "planck_comparison_path": planck_comp_path
    }


def process_all_existing_checkpoints(checkpoints_dir: str = "checkpoints"):
    """Scans and processes all existing eon_*.npz checkpoints, generating full observational plots."""
    import glob
    npz_files = sorted(glob.glob(os.path.join(checkpoints_dir, "eon_*.npz")))
    if not npz_files:
        print(f"No eon_*.npz files found in '{checkpoints_dir}'.")
        return

    print(f"• Found {len(npz_files)} historical eon checkpoints in '{checkpoints_dir}'. Processing observational reports...")
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
            engine.scale_factor = float(data['scale_factor'])
            engine.eon = int(data['eon'])
            engine.total_steps = int(data.get('total_steps', 0))
            if 'tau_eon_start' in data:
                engine.tau_eon_start = data['tau_eon_start']

            print(f"\n--- Processing {os.path.basename(npz_path)} (Eon {engine.eon}) ---")
            metrics = generate_eon_observational_report(engine, output_dir=snapshots_dir)
            print(f"  * C2/C3 Ratio: {metrics['ratio_C2_C3']} | Predicted H0: {metrics['h0_predicted']} km/s/Mpc | Chi2: {metrics['planck_chi2']}")
        except Exception as ex:
            print(f"  * Error processing {npz_path}: {ex}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Reotransductor vs. ESA Planck 2018 Observational Validation")
    parser.add_argument("--checkpoint-dir", type=str, default="checkpoints", help="Directory with latest checkpoint or snapshot")
    parser.add_argument("--steps", type=int, default=0, help="Additional integration steps to simulate before evaluation")
    parser.add_argument("--from-scratch", action="store_true", help="Start from blank initial conditions rather than active checkpoint")
    parser.add_argument("--process-all", action="store_true", help="Process and generate plots for all historical eon_*.npz checkpoints")
    parser.add_argument("--output", type=str, default="assets/planck_comparison_spectrum.png", help="Output PNG path")
    args = parser.parse_args()

    if args.process_all:
        process_all_existing_checkpoints(checkpoints_dir=args.checkpoint_dir)
    else:
        run_planck_comparison(
            checkpoint_path=args.checkpoint_dir,
            auto_resume=not args.from_scratch,
            steps=args.steps,
            output_fig=args.output
        )
