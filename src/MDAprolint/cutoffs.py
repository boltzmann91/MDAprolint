import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

# Import your optimized distance calculator
from .distances import get_lipid_protein_distances

def calculate_optimal_cutoffs(u, prot_sel, lip_sel, step=1, max_dist=15.0, bins=200, sigma=3.0):
    """
    Calculates distances, fits a smooth curve, and identifies the first peak 
    and first valley to dictate the lower and upper dual-cutoffs.
    
    Returns:
    - tuple: (lower_cutoff, upper_cutoff) in Angstroms.
    """
    
    # Automated Directory and Filename Handling
    output_dir = "cutoffs"
    os.makedirs(output_dir, exist_ok=True)
    formatted_lip_name = lip_sel.replace(" ", "_")
    save_path = os.path.join(output_dir, f"{formatted_lip_name}_dual_cutoffs.png")

    print("--- Step 1: Calculating distances for cutoff analysis ---")
    lip_resids, dist_array, nearest_resids = get_lipid_protein_distances(u, prot_sel, lip_sel, step=step)
    
    print("--- Step 2: Fitting curve to find optimal dual cutoffs ---")
    all_distances = dist_array.flatten()
    counts, bin_edges = np.histogram(all_distances, bins=bins, range=(0, max_dist), density=True)
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    smoothed_curve = gaussian_filter1d(counts, sigma=sigma)
    
    # Find Peaks and Valleys
    peaks, _ = find_peaks(smoothed_curve, prominence=0.005)
    valleys, _ = find_peaks(-smoothed_curve, prominence=0.001)
    
    if len(peaks) == 0 or len(valleys) == 0:
        raise ValueError("Could not detect clear peaks/valleys.")
        
    # --- DUAL CUTOFF LOGIC ---
    # Lower Cutoff = First Main Peak
    first_main_peak_idx = peaks[0]
    lower_cutoff = bin_centers[first_main_peak_idx]
    
    # Upper Cutoff = First Valley after the Main Peak
    valid_valleys = [v for v in valleys if v > first_main_peak_idx]
    if len(valid_valleys) == 0:
        raise ValueError("Could not find a valley after the first peak.")
        
    first_valley_idx = valid_valleys[0]
    upper_cutoff = bin_centers[first_valley_idx]
    
    # --- Create the Plot ---
    fig, ax = plt.subplots(figsize=(8, 5))
    
    ax.hist(all_distances, bins=bins, range=(0, max_dist), density=True, 
            color='lightgray', edgecolor='none', label='Raw Distribution')
    ax.plot(bin_centers, smoothed_curve, color='blue', linewidth=2.5, label='Fitted Curve')
    
    # Draw Dual Cutoff Lines
    ax.axvline(lower_cutoff, color='green', linestyle='--', linewidth=2,
               label=f'Lower Cutoff (Peak): {lower_cutoff:.2f} Å')
    ax.axvline(upper_cutoff, color='red', linestyle='--', linewidth=2,
               label=f'Upper Cutoff (Valley): {upper_cutoff:.2f} Å')
    
    # Shade the "Rattling" Buffer Zone
    ax.axvspan(lower_cutoff, upper_cutoff, color='yellow', alpha=0.2, label='Buffer Zone')
    
    ax.set_title("Automated Dual-Cutoff Detection")
    ax.set_xlabel("Minimum Distance to Protein (Å)")
    ax.set_ylabel("Probability Density")
    ax.set_xlim(0, max_dist)
    ax.legend(loc='upper right')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close(fig)
    
    print(f"--> Lower Cutoff (Interaction Starts): {lower_cutoff:.2f} Å")
    print(f"--> Upper Cutoff (Interaction Ends):   {upper_cutoff:.2f} Å")
    print(f"--> Diagnostic plot saved to: {save_path}\n")
    
    return lower_cutoff, upper_cutoff