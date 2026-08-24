import os
import sys
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from observational.galaxy_simulator import run_single_galaxy_simulation

def main():
    sparc_path = "data/sparc_2020/sparc_rotation_curves.json"
    if not os.path.exists(sparc_path):
        print(f"File not found: {sparc_path}")
        return
        
    with open(sparc_path, 'r') as f:
        sparc_raw = json.load(f)
        
    galaxies = list(sparc_raw.get("galaxies", {}).keys())
    total_gals = len(galaxies)
    
    out_dir = "assets/galaxies"
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"Starting batch analysis for {total_gals} galaxies...")
    
    results_list = []
    reot_wins = 0
    nfw_wins = 0
    
    t0 = time.time()
    
    # Hide output of the single runs using a simple redirect or just let it print
    # Let's let it print to keep the log active
    for i, gal in enumerate(galaxies):
        print(f"[{i+1}/{total_gals}] Processing {gal}...")
        try:
            fig_path = os.path.join(out_dir, f"galaxy_{gal.lower()}_rot_curve.png")
            res = run_single_galaxy_simulation(galaxy_name=gal, grid_size=64, output_fig=fig_path)
            results_list.append(res)
            if "Reotransductor" in res["preferred_model"]:
                reot_wins += 1
            else:
                nfw_wins += 1
        except Exception as e:
            print(f"  Error processing {gal}: {e}")
            
    t1 = time.time()
    print(f"\n=======================================================")
    print(f"Batch completed in {t1 - t0:.2f} seconds.")
    
    # Generate Markdown Report
    md_path = "assets/sparc_population_results.md"
    
    total_processed = len(results_list)
    reot_pct = (reot_wins / total_processed * 100) if total_processed > 0 else 0
    nfw_pct = (nfw_wins / total_processed * 100) if total_processed > 0 else 0
    
    with open(md_path, 'w') as f:
        f.write("# SPARC Population Analysis (175 Galaxies)\n\n")
        f.write("## Overall Statistical Validation\n")
        f.write(f"- **Total Galaxies Analyzed:** {total_processed}\n")
        f.write(f"- **Reotransductor Cored Profile Preferred:** {reot_wins} galaxies ({reot_pct:.1f}%)\n")
        f.write(f"- **NFW Cuspy Model Preferred:** {nfw_wins} galaxies ({nfw_pct:.1f}%)\n\n")
        
        f.write("## Galaxy Metrics Data\n")
        f.write("| Galaxy | Type | Dist (Mpc) | Reotransductor $\\chi^2_\\nu$ | NFW $\\chi^2_\\nu$ | Preferred Model |\n")
        f.write("|--------|------|------------|-----------------------------|-----------------|-----------------|\n")
        
        for r in results_list:
            reot_chi2 = r.get("red_chi2_reotransductor", "N/A")
            nfw_chi2 = r.get("red_chi2_nfw", "N/A")
            pref = r.get("preferred_model", "N/A")
            name = r.get("galaxy_name", "Unknown")
            gtype = r.get("galaxy_type", "Unknown")
            dist = r.get("distance_mpc", "Unknown")
            f.write(f"| {name} | {gtype} | {dist} | {reot_chi2} | {nfw_chi2} | {pref} |\n")
            
    print(f"Markdown report generated at {md_path}")
    print(f"=======================================================\n")
    
if __name__ == '__main__':
    main()
