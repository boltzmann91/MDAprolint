import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import distance_array
import numpy as np

class LipidProteinDistanceAnalysis(AnalysisBase):
    """
    Calculates the minimum distance from each lipid atom to the protein per frame,
    and identifies the nearest protein residue.
    """
    def __init__(self, universe, prot_sel, lip_sel, **kwargs):
        super().__init__(universe.trajectory, **kwargs)
        self.prot_atoms = universe.select_atoms(prot_sel)
        self.lip_atoms = universe.select_atoms(lip_sel)
        
        # Pre-extract resids for mapping
        self.lip_resids = self.lip_atoms.resids
        self.prot_resids = self.prot_atoms.resids

    def _prepare(self):
        # We will track two separate lists of 1D arrays
        self.all_distances = []
        self.all_nearest_resids = []

    def _single_frame(self):
        # Calculate the distance matrix (n_protein_atoms, n_lipid_atoms)
        dist_matrix = distance_array(
            self.prot_atoms.positions,
            self.lip_atoms.positions,
            box=self._ts.dimensions,
            backend='OpenMP'
        )
        
        # 1. Get the shortest distance for each lipid atom
        min_dists_per_lipid = dist_matrix.min(axis=0)
        
        # 2. Get the index of the protein atom that is closest
        min_indices = dist_matrix.argmin(axis=0)
        
        # 3. Map that index back to the actual protein residue ID
        nearest_prot_resids = self.prot_resids[min_indices]
        
        # Store both for this frame
        self.all_distances.append(min_dists_per_lipid)
        self.all_nearest_resids.append(nearest_prot_resids)

    def _conclude(self):
        # Convert both lists to 2D numpy arrays
        # Shapes: (n_frames, n_lipid_atoms)
        self.results_distances = np.array(self.all_distances)
        self.results_nearest_resids = np.array(self.all_nearest_resids)


def get_lipid_protein_distances(u, prot_sel, lip_sel, step=1):
    """
    Wrapper function to execute continuous distance and nearest-residue calculation.
    
    Returns:
    - lipid_resids (1D numpy array): The specific lipid resids mapping to the columns.
    - dist_array (2D numpy array): Distances over time.
    - nearest_resid_array (2D numpy array): Closest protein resids over time.
    """
    print(f"Calculating continuous distances over {u.trajectory.n_frames} frames...")
    analysis = LipidProteinDistanceAnalysis(u, prot_sel, lip_sel)
    analysis.run(step=step, verbose=True)
    
    lipid_resids = analysis.lip_atoms.resids
    
    return lipid_resids, analysis.results_distances, analysis.results_nearest_resids