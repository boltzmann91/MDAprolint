import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import distance_array
import numpy as np

class LipidProteinDistanceAnalysis(AnalysisBase):
    """
    Calculates the minimum distance from each lipid atom to the protein per frame.
    """
    def __init__(self, universe, prot_sel, lip_sel, **kwargs):
        super().__init__(universe.trajectory, **kwargs)
        self.prot_atoms = universe.select_atoms(prot_sel)
        self.lip_atoms = universe.select_atoms(lip_sel)

    def _prepare(self):
        self.all_distances = []

    def _single_frame(self):
        dist_matrix = distance_array(
            self.prot_atoms.positions,
            self.lip_atoms.positions,
            box=self._ts.dimensions,
            backend='OpenMP'
        )
        
        min_dists_per_lipid = dist_matrix.min(axis=0)
        self.all_distances.append(min_dists_per_lipid)

    def _conclude(self):
        self.results = np.array(self.all_distances)

def get_lipid_protein_distances(u, prot_sel, lip_sel, step=1):
    """
    Wrapper function to execute the continuous distance calculation.
    """
    print(f"Calculating continuous distances over {u.trajectory.n_frames} frames...")
    analysis = LipidProteinDistanceAnalysis(u, prot_sel, lip_sel)
    analysis.run(step=step, verbose=True)
    
    lipid_resids = analysis.lip_atoms.resids
    
    return analysis.results, lipid_resids
