import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

from .distances import get_lipid_protein_distances

def calculate_optimal_cutoffs(u, prot_sel, lip_sel, step=1, max_dist=15.0, bins=200, sigma=1.0):
    """
    Calculates distances, generates a density curve, and uses spline calculus
    to detect the optimal dual-cutoffs. 
    
    Parameters:
    - u: MDAnalysis Universe object.
    - prot_sel (str): Selection string for the protein.
    - lip_sel (str): Selection string for the lipids.
    - step (int): Frame step size for reading the trajectory.
    - max_dist (float): Maximum distance to analyze on the x-axis.
    - bins (int): Number of bins for the high-res histogram.
    - smoothing (float): Smoothing factor for the spline fit.
    - save_plot (bool): If True, generates and saves a diagnostic plot (Default: False).
    
    Returns:
    - tuple: (lower_cutoff, upper_cutoff) in Angstroms.
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    
    print("--- Step 1: Calculating distances for cutoff analysis ---")
    lip_resids, dist_array, nearest_resids = get_lipid_protein_distances(u, prot_sel, lip_sel, step=step)
    
    print("--- Step 2: Applying spline calculus to find exact cutoffs ---")
    all_distances = dist_array.flatten()
    
    # Generate the raw numerical inputs (X and Y)
    counts, bin_edges = np.histogram(all_distances, bins=bins, range=(0, max_dist), density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    smoothed_curve = gaussian_filter1d(counts, sigma=sigma)
    
    # Find Peaks and Valleys
    peaks, _ = find_peaks(smoothed_curve, prominence=0.005)
    valleys, _ = find_peaks(-smoothed_curve, prominence=0.0005)
    
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
                color='lightgray', edgecolor='none', label='Raw Histogram')
    ax.plot(x_smooth, y_smooth, color='blue', linewidth=2.5, label='Fitted Spline')
        
    ax.axvline(lower_cutoff, color='green', linestyle='--', linewidth=2,
                   label=f'Lower Cutoff (Peak): {lower_cutoff:.2f} Å')
    ax.axvline(upper_cutoff, color='red', linestyle='--', linewidth=2,
                   label=f'Upper Cutoff (Valley): {upper_cutoff:.2f} Å')
        
    ax.axvspan(lower_cutoff, upper_cutoff, color='yellow', alpha=0.2, label='Buffer Zone')
        
    ax.set_title("Spline Derivative Dual-Cutoff Detection")
    ax.set_xlabel("Minimum Distance to Protein (Å)")
    ax.set_ylabel("Probability Density")
    ax.set_xlim(0, max_dist)
    ax.legend(loc='upper right')
    ax.grid(axis='y', alpha=0.3)
        
    plt.tight_layout()
    plt.show()
    if save_plot:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"--> Diagnostic plot saved to: {save_path}\n")
    else:
        print("\n") # Add a clean newline if not printing the save path
        plt.close(fig)
    return lower_cutoff, upper_cutoff