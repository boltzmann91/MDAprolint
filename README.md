***

# MDAprolint 
**MDAnalysis to examine Protein-Lipid Interactions**

`MDAprolint` is a high-performance Python package designed to analyze lipid-protein interactions in Molecular Dynamics (MD) simulations. Built on top of `MDAnalysis` and principles of Prolint, `MDAprolint` is accelerated by OpenMP and provides an end-to-end pipeline for calculating lipid residence times, filtering out MD noise, and visualizing 3D binding hotspots.

## ✨ Key Features
* **Automated Dual-Cutoff Detection:** Empirically calculates primary annular shell boundaries from your trajectory data using distance distributions and Gaussian smoothing.
* **Rigorous Noise Filtering:** Uses a dual-cutoff state machine to eliminate false binding events caused by thermal "rattling" and second-shell "drive-bys."
* **OpenMP Acceleration:** Instantly calculates massive distance matrices for both Atomic and Center of Mass (COM) selections.
* **Kinetics & Lifetimes:** Outputs clean Pandas DataFrames of continuous binding durations and docking sites.
* **3D Hotspot Visualization:** Automatically maps per-residue occupancy data to the B-factor column of a PDB for instant visualization in PyMOL or VMD.

---

## ⚙️ Installation

We recommend installing `MDAprolint` inside an isolated virtual environment (e.g., using `conda` or `venv`).

**Preferred Method (Stable Release)**  
Install the latest stable version directly from PyPI:
```bash
pip install MDAprolint
```

**Development Version**  
If you want the absolute latest features or wish to modify the code yourself, you can install the development version directly from the repository:
```bash
git clone https://github.com/boltzmann91/MDAprolint.git
cd MDAprolint
pip install -e .
```
*(Note: Both methods automatically install all required dependencies, including MDAnalysis, Pandas, NumPy, Matplotlib, and SciPy).*

---

## 🚀 Quick Start Workflow

Here is a suggested workflow that takes you from a raw trajectory to a 3D PyMOL hotspot structure in just a few lines of code.

```python
import MDAnalysis as mda
from MDAprolint import (
    calculate_optimal_cutoffs,
    get_interacting_lipids,
    calculate_residence_times,
    calculate_occupancy,
    plot_occupancy,
    export_occupancy_bfactors
)

# 1. Load your trajectory
u = mda.Universe("system.psf", "trajectory.dcd")

# Define selections (Use 'not type H*' for best AA-MD accuracy)
prot_sel = "protein and not type H*"
lip_sel = "resname POPC and not type H*"
lipid_name = "POPC"

# 2. Automate the Physics (Calculate exact cutoffs from the data)
# save_plot=True generates a diagnostic distribution graph in the 'cutoffs/' folder.
lower_cut, upper_cut = calculate_optimal_cutoffs(
    u, prot_sel, lip_sel, 
    step=10, 
    save_plot=True
)

# 3. Run the Dual-Cutoff Interaction Engine
interactions = get_interacting_lipids(
    u, prot_sel, lip_sel, 
    lower_cutoff=lower_cut, 
    upper_cutoff=upper_cut,
    step=1 
)

# 4. Extract Kinetics (Residence Times)
# Assumes 10 ps between saved trajectory frames
df_residence = calculate_residence_times(interactions, frame_time=10.0)

# Filter out fast collisions (< 500 ps) to find stable binders
df_stable = df_residence[df_residence['duration'] >= 500.0]

# 5. Calculate & Visualize Occupancy
df_occ = calculate_occupancy(interactions)
plot_occupancy(df_occ) # Saves a bar chart to the current directory

# 6. Export 3D PyMOL Hotspots
# Passing a clean (PSF, PDB) prevents periodic boundary wrapping issues.
export_occupancy_bfactors(
    input_struct=("clean.psf", "clean.pdb"), 
    prot_sel=prot_sel, 
    occupancy_df=df_occ, 
    lipid_name=lipid_name
)
```

---

## 🎨 Visualizing Hotspots in PyMOL

Running the workflow above will automatically generate a new folder called `occupancy-PDBs/` containing your customized `.pdb` file.

To visualize the lipid binding sites:
1. Open the `.pdb` file in **PyMOL**.
2. In the PyMOL command line, type the following commands:
```pymol
hide all
show surface
spectrum b, white_red
```
Areas colored **white** indicate 0% lipid occupancy, while areas colored **dark red** indicate high-affinity binding pockets.

---

## 🧬 Advanced Usage: Rigid Sterols (COM tracking)

For rigid ring systems like **Cholesterol** or **Ergosterol**, measuring interactions via Center of Mass (COM) is often preferred to bypass headgroup tilting artifacts.

`MDAprolint` supports this natively. Simply add `use_com=True` to the relevant functions and select the ring atoms:

```python
# Select only the rigid sterol rings
sterol_rings = "resname ERG and name C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17"

# 1. Calculate Cutoffs using COM
lower_cut, upper_cut = calculate_optimal_cutoffs(
    u, prot_sel, sterol_rings, 
    use_com=True, 
    save_plot=True
)

# 2. Track interactions using COM
interactions = get_interacting_lipids(
    u, prot_sel, sterol_rings, 
    lower_cutoff=lower_cut, upper_cutoff=upper_cut, 
    use_com=True
)
```
