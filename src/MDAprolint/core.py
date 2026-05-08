import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import distance_array
import pandas as pd
import numpy as np

# =============================================================================
# Function 1: Find Interacting Lipids (AnalysisBase)
# =============================================================================
class InteractingLipidsAnalysis(AnalysisBase):
    def __init__(self, universe, prot_sel, lip_sel, cutoff, **kwargs):
        super().__init__(universe.trajectory, **kwargs)
        self.prot_atoms = universe.select_atoms(prot_sel)
        self.lip_atoms = universe.select_atoms(lip_sel)
        self.cutoff = cutoff
        self.lip_resids = self.lip_atoms.resids

    def _prepare(self):
        self.interactions_per_frame = []

    def _single_frame(self):
        dist_matrix = distance_array(
            self.prot_atoms.positions,
            self.lip_atoms.positions,
            box=self._ts.dimensions,
            backend='OpenMP'
        )
        
        interacting_mask = (dist_matrix <= self.cutoff).any(axis=0)
        interacting_resids = set(self.lip_resids[interacting_mask])
        
        self.interactions_per_frame.append({
            'frame': self._ts.frame,
            'lipids': interacting_resids
        })

    def _conclude(self):
        self.results = self.interactions_per_frame

def get_interacting_lipids(u, prot_sel, lip_sel, cutoff, step=1):
    print(f"Finding interactions over {u.trajectory.n_frames} frames...")
    analysis = InteractingLipidsAnalysis(u, prot_sel, lip_sel, cutoff)
    analysis.run(step=step, verbose=True)
    return analysis.results


# =============================================================================
# Function 2: Calculate Residence Times
# =============================================================================
def calculate_residence_times(interactions_per_frame, frame_time):
    active_interactions = {} 
    all_events = []          
    
    print("Calculating residence and cumulative times...")
    for step_index, frame_data in enumerate(interactions_per_frame):
        current_lipids = frame_data['lipids']
        current_time = step_index * frame_time
        
        for lipid in list(active_interactions.keys()):
            if lipid not in current_lipids:
                start_time = active_interactions.pop(lipid)
                end_time = current_time
                duration = end_time - start_time
                
                all_events.append({
                    'lipid_id': lipid, 'start_time': start_time,
                    'end_time': end_time, 'duration': duration
                })
                
        for lipid in current_lipids:
            if lipid not in active_interactions:
                active_interactions[lipid] = current_time
                
    final_time = (len(interactions_per_frame) - 1) * frame_time
    for lipid, start_time in active_interactions.items():
        all_events.append({
            'lipid_id': lipid, 'start_time': start_time,
            'end_time': final_time, 'duration': final_time - start_time
        })
        
    df = pd.DataFrame(all_events)
    if not df.empty:
        df = df.sort_values(by=['lipid_id', 'start_time']).reset_index(drop=True)
    return df


# =============================================================================
# NEW Function 3: Calculate Continuous Distances (AnalysisBase)
# =============================================================================
class LipidProteinDistanceAnalysis(AnalysisBase):
    """
    Calculates the minimum distance from each lipid atom to the protein per frame.
    """
    def __init__(self, universe, prot_sel, lip_sel, **kwargs):
        super().__init__(universe.trajectory, **kwargs)
        self.prot_atoms = universe.select_atoms(prot_sel)
        self.lip_atoms = universe.select_atoms(lip_sel)

    def _prepare(self):
        # We will store a list of 1D numpy arrays, one for each frame
        self.all_distances = []

    def _single_frame(self):
        # dist_matrix shape: (n_protein_atoms, n_lipid_atoms)
        dist_matrix = distance_array(
            self.prot_atoms.positions,
            self.lip_atoms.positions,
            box=self._ts.dimensions,
            backend='OpenMP'
        )
        
        # .min(axis=0) finds the shortest distance down the protein column.
        # Resulting shape: (n_lipid_atoms,) - i.e. one distance per lipid atom.
        min_dists_per_lipid = dist_matrix.min(axis=0)
        
        self.all_distances.append(min_dists_per_lipid)

    def _conclude(self):
        # Convert the list of 1D arrays into a 2D numpy array
        # Final shape: (n_frames, n_lipid_atoms)
        self.results = np.array(self.all_distances)


def get_lipid_protein_distances(u, prot_sel, lip_sel, step=1):
    """
    Wrapper function to execute the continuous distance calculation.
    
    Returns:
    - 2D Numpy Array of shape (frames, lipid_atoms)
    - 1D Numpy Array mapping the columns to specific lipid resids
    """
    print(f"Calculating continuous distances over {u.trajectory.n_frames} frames...")
    analysis = LipidProteinDistanceAnalysis(u, prot_sel, lip_sel)
    analysis.run(step=step, verbose=True)
    
    # We also return the resids so you know which column belongs to which lipid!
    lipid_resids = analysis.lip_atoms.resids
    
    return analysis.results, lipid_resids


# =============================================================================
# Main Execution Block
# =============================================================================
if __name__ == "__main__":
    
    topology_file = "system.pdb"
    trajectory_file = "trajectory.xtc"
    
    # u = mda.Universe(topology_file, trajectory_file)
    
    protein_selection = "protein and not type H*"
    lipid_selection = "resname POPC and name P"
    cutoff_distance = 5.0 
    frame_time_ps = 10.0 
    
    # --- 1. Boolean Interaction Search & Residence Times ---
    # interactions = get_interacting_lipids(u, protein_selection, lipid_selection, cutoff_distance)
    # df_res_times = calculate_residence_times(interactions, frame_time_ps)
    # print(df_res_times.head())
    
    # --- 2. NEW: Continuous Distance Tracking ---
    # dist_array, lip_resids = get_lipid_protein_distances(u, protein_selection, lipid_selection)
    
    # print(f"Distance Array Shape: {dist_array.shape}")
    # Example Output: Distance Array Shape: (1001, 250) -> 1001 frames, 250 lipid atoms
    
    # Example: How to find the distance of the 5th lipid over time
    # fifth_lipid_resid = lip_resids[4]
    # fifth_lipid_distances = dist_array[:, 4]
    # print(f"Distances for Lipid {fifth_lipid_resid}: {fifth_lipid_distances}")


