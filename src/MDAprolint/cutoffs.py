import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

from .distances import get_lipid_protein_distances

def find_cutoffs_from_curve(x, y, smoothing_factor=0.01):
    """
    Takes raw numerical arrays (X, Y), applies a continuous spline fit, 
    and uses calculus to locate the FIRST significant peak and its subsequent valley.
    """
    from scipy.interpolate import UnivariateSpline
    import numpy as np
    
    # 1. Fit a continuous mathematical spline
    spline = UnivariateSpline(x, y, s=smoothing_factor)
    
    # 2. Calculate the derivative functions
    dy = spline.derivative(n=1)  
    ddy = spline.derivative(n=2) 
    
    # 3. Create an ultra-dense grid to evaluate zero-crossings
    x_dense = np.linspace(x.min(), x.max(), 10000)
    dy_dense = dy(x_dense)
    zero_crossings = np.where(np.diff(np.sign(dy_dense)))[0]
    
    if len(zero_crossings) == 0:
        raise ValueError("Could not find any peaks or valleys.")
        
    roots = x_dense[zero_crossings]
    
    # 4. Sort the roots into Peaks and Valleys
    peaks = [r for r in roots if ddy(r) < 0]
    valleys = [r for r in roots if ddy(r) > 0]
    
    if not peaks:
        raise ValueError("Could not find a valid peak.")
        
    # --- THE UPDATED LOGIC FOR TALLER SECOND PEAKS ---
    
    # Calculate the height (density) of every peak
    peak_densities = np.array([spline(p) for p in peaks])
    max_density = np.max(peak_densities)
    
    # Filter: A peak must be at least 15% as tall as the highest peak to be considered "real"
    # This ignores tiny thermal ripples at the far left of the graph
    threshold = 0.15 * max_density
    
    valid_peaks = [p for p, d in zip(peaks, peak_densities) if d >= threshold]
    
    if not valid_peaks:
        raise ValueError("No significant peaks found above the noise threshold.")
    
    # The Lower Cutoff is the FIRST valid peak (the one physically closest to the protein)
    lower_cutoff = valid_peaks[0]
    
    # The Upper Cutoff is the FIRST valley AFTER this lower cutoff
    valid_valleys = [v for v in valleys if v > lower_cutoff]
    
    if not valid_valleys:
        raise ValueError("Could not find a valid valley after the first peak.")
        
    upper_cutoff = valid_valleys[0]
    
    return lower_cutoff, upper_cutoff, spline


def calculate_optimal_cutoffs(u, prot_sel, lip_sel, step=1, max_dist=15.0, bins=200, smoothing=0.01, save_plot=False):
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
    
    # Use the decoupled math function to get the exact cutoffs!
    lower_cutoff, upper_cutoff, fitted_spline = find_cutoffs_from_curve(
        x=bin_centers, 
        y=counts, 
        smoothing_factor=smoothing
    )
    
    print(f"--> Lower Cutoff: {lower_cutoff:.2f} Å")
    print(f"--> Upper Cutoff: {upper_cutoff:.2f} Å")
    
    # --- OPTIONAL PLOTTING BLOCK ---
    
    output_dir = "cutoffs"
    os.makedirs(output_dir, exist_ok=True)
    formatted_lip_name = lip_sel.replace(" ", "_").replace("*", "")
    save_path = os.path.join(output_dir, f"{formatted_lip_name}_dual_cutoffs.png")
        
    # Generate an ultra-smooth line for plotting
    x_smooth = np.linspace(0, max_dist, 1000)
    y_smooth = fitted_spline(x_smooth)
        
    # Create the Plot
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