import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

# Import your optimized distance calculator from the sibling module
from .distances import get_lipid_protein_distances

def calculate_optimal_cutoff(u, prot_sel, lip_sel, step=1, save_path="cutoff_analysis.png", max_dist=15.0, bins=200, sigma=3.0):
    """
    Calculates distances, fits a smooth curve, and identifies the first 
    annular shell boundary (valley) to suggest an interaction cutoff.
    
    Parameters:
    - u: MDAnalysis Universe object.
    - prot_sel (str): Selection string for the protein.
    - lip_sel (str): Selection string for the lipids.
    - step (int): Frame step size for reading the trajectory.
    - save_path (str): Filename to save the diagnostic plot.
    - max_dist (float): Maximum distance to analyze on the x-axis.
    - bins (int): Number of bins for the high-res histogram.
    - sigma (float): Smoothing factor for the Gaussian curve fit.
    
    Returns:
    - float: The suggested cutoff distance in Angstroms.
    """
    print("--- Step 1: Calculating distances for cutoff analysis ---")
    # Run the distance calculation under the hood
    lip_resids, dist_array, nearest_resids = get_lipid_protein_distances(u, prot_sel, lip_sel, step=step)
    
    print("--- Step 2: Fitting curve to find optimal cutoff ---")
    
    # 1. Flatten data and calculate high-resolution histogram
    all_distances = dist_array.flatten()
    counts, bin_edges = np.histogram(all_distances, bins=bins, range=(0, max_dist), density=True)
    
    # Get the center x-coordinates of each bin
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # 2. Fit a smooth curve using a Gaussian filter
    smoothed_curve = gaussian_filter1d(counts, sigma=sigma)
    
    # 3. Find Peaks and Valleys
    peaks, _ = find_peaks(smoothed_curve, prominence=0.005)
    valleys, _ = find_peaks(-smoothed_curve, prominence=0.001)
    
    if len(peaks) == 0 or len(valleys) == 0:
        raise ValueError("Could not detect clear peaks/valleys. Adjust smoothing 'sigma' or check your trajectory.")
        
    # We want the FIRST valley that occurs AFTER the FIRST main peak
    first_main_peak_idx = peaks[0]
    valid_valleys = [v for v in valleys if v > first_main_peak_idx]
    
    if len(valid_valleys) == 0:
        raise ValueError("Could not find a valley after the first peak.")
        
    first_valley_idx = valid_valleys[0]
    suggested_cutoff = bin_centers[first_valley_idx]
    
    # 4. Create the Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    
    ax.hist(all_distances, bins=bins, range=(0, max_dist), density=True, 
            color='lightgray', edgecolor='none', label='Raw Distribution')
    
    ax.plot(bin_centers, smoothed_curve, color='blue', linewidth=2.5, label='Fitted Curve')
    
    ax.plot(bin_centers[first_main_peak_idx], smoothed_curve[first_main_peak_idx], 'go', label="1st Shell Peak")
    ax.plot(bin_centers[first_valley_idx], smoothed_curve[first_valley_idx], 'ro', label="1st Shell Boundary")
    
    ax.axvline(suggested_cutoff, color='red', linestyle='--', linewidth=2,
               label=f'Suggested Cutoff: {suggested_cutoff:.2f} Å')
    
    ax.set_title("Automated Annular Shell Boundary Detection")
    ax.set_xlabel("Minimum Distance to Protein (Å)")
    ax.set_ylabel("Probability Density")
    ax.set_xlim(0, max_dist)
    ax.legend(loc='upper right')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    print(f"--> Suggested Cutoff successfully calculated: {suggested_cutoff:.2f} Å")
    print(f"--> Diagnostic plot saved to: {save_path}\n")
    
    return suggested_cutoff