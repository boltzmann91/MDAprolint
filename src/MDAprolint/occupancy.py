import os
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda

def calculate_occupancy(interactions_per_frame):
    """
    Calculates the per-residue lipid occupancy percentage.
    
    Parameters:
    - interactions_per_frame: List of dicts (output from get_interacting_lipids).
    
    Returns:
    - Pandas DataFrame with protein_resid and occupancy_percent.
    """
    print("Calculating per-residue occupancy...")
    total_frames = len(interactions_per_frame)
    residue_counts = {}

    for frame_data in interactions_per_frame:
        # frame_data['lipids'] is a dict of {lipid_resid: protein_resid}
        # .values() gives us all the protein residues bound in this specific frame
        occupied_resids = set(frame_data['lipids'].values())
        
        for res in occupied_resids:
            residue_counts[res] = residue_counts.get(res, 0) + 1

    # Convert the raw counts into a clean DataFrame
    data = []
    for res, count in residue_counts.items():
        data.append({
            'protein_resid': res,
            'occupancy_frames': count,
            'occupancy_percent': (count / total_frames) * 100.0
        })

    df = pd.DataFrame(data)
    if not df.empty:
        # Sort by residue number so it's in sequential order along the protein backbone
        df = df.sort_values(by='protein_resid').reset_index(drop=True)
        
    return df


def plot_occupancy(occupancy_df, save_path="occupancy_plot.png"):
    """
    Generates a bar chart of lipid occupancy along the protein sequence.
    """
    fig, ax = plt.subplots(figsize=(12, 4))
    
    # Plot the bar chart
    ax.bar(occupancy_df['protein_resid'], occupancy_df['occupancy_percent'], color='salmon', width=1.0)
    
    ax.set_title("Per-Residue Lipid Occupancy")
    ax.set_xlabel("Protein Residue ID")
    ax.set_ylabel("Occupancy (%)")
    ax.set_ylim(0, 100) # Occupancy goes from 0 to 100%
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"--> Occupancy plot saved to: {save_path}")

import os
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda

# ... (calculate_occupancy and plot_occupancy remain unchanged) ...

def export_occupancy_bfactors(input_struct, prot_sel, occupancy_df, lipid_name="LIPID", save_path=None, frame=0):
    """
    Exports a PDB where the B-factor column is replaced by Occupancy %.
    Automatically saves to an 'occupancy-PDBs/' directory.
    """
    
    # Automated Directory and Filename Handling
    if save_path is None:
        output_dir = "occupancy-PDBs"
        os.makedirs(output_dir, exist_ok=True)
        formatted_prot = prot_sel.replace(" ", "_").replace("*", "")
        filename = f"{formatted_prot}-{lipid_name}-hotspots.pdb"
        save_path = os.path.join(output_dir, filename)

    #  Intelligently handle the input structure
    if isinstance(input_struct, mda.Universe):
        print(f"Using provided Universe (rewinding to frame {frame})...")
        u = input_struct
        u.trajectory[frame]
    elif isinstance(input_struct, str):
        print(f"Loading reference structure from file: {input_struct}...")
        u = mda.Universe(input_struct)
    elif isinstance(input_struct, (list, tuple)):
        print(f"Loading reference structure from files: {input_struct}...")
        u = mda.Universe(*input_struct)
    else:
        raise ValueError("input_struct must be an mda.Universe, a filename string, or a tuple of filenames.")
    
    # Check if the entire Universe has the 'tempfactors' (B-factors) attribute.
    # If not, add it to the universe. MDAnalysis initializes them all to 0.0 automatically.
    if not hasattr(u.atoms, 'tempfactors'):
        print("B-factors not found in input topology. Initializing them to 0.0...")
        u.add_TopologyAttr('tempfactors')
    
    
    # (Just in case the loaded PDB had old, leftover experimental B-factors)
    prot = u.select_atoms(prot_sel)
    prot.tempfactors = 0.0
    

    for _, row in occupancy_df.iterrows():
        res_id = int(row['protein_resid'])
        occ = row['occupancy_percent']
        
        # Select the atoms of this specific residue inside our protein group
        res_atoms = prot.select_atoms(f"resid {res_id}")
        
        # Only assign if the residue exists in the selection to prevent errors
        if len(res_atoms) > 0:
            res_atoms.tempfactors = occ

    prot.write(save_path)
    print(f"--> PyMOL PDB exported to: {save_path}")
    print("    (In PyMOL, type: 'spectrum b, white_red' to see the hotspots!)")