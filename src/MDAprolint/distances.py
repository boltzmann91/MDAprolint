import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import distance_array
import numpy as np

class LipidProteinDistanceAnalysis(AnalysisBase):
    def __init__(self, universe, prot_sel, lip_sel, use_com=False, **kwargs):
        super().__init__(universe.trajectory, **kwargs)
        self.prot_atoms = universe.select_atoms(prot_sel)
        self.lip_atoms = universe.select_atoms(lip_sel)
        self.use_com = use_com
        
        # Pre-extract resids for mapping
        if self.use_com:
            # If using COM, we get 1 coordinate per lipid molecule
            self.lip_resids = self.lip_atoms.residues.resids
        else:
            # If using atoms, we get 1 coordinate per lipid atom
            self.lip_resids = self.lip_atoms.resids
            
        self.prot_resids = self.prot_atoms.resids

    def _prepare(self):
        self.all_distances = []
        self.all_nearest_resids = []

    def _single_frame(self):
        # 1. Get Lipid Positions (Atoms vs. Center of Mass)
        if self.use_com:
            # Calculates COM for the selected atoms, grouped by residue
            lip_positions = self.lip_atoms.center_of_mass(compound='residues')
        else:
            lip_positions = self.lip_atoms.positions

        # 2. Calculate Distance Matrix
        dist_matrix = distance_array(
            self.prot_atoms.positions,
            lip_positions,
            box=self._ts.dimensions,
            backend='OpenMP'
        )
        
        min_dists = dist_matrix.min(axis=0)
        min_indices = dist_matrix.argmin(axis=0)
        nearest_prot = self.prot_resids[min_indices]
        
        self.all_distances.append(min_dists)
        self.all_nearest_resids.append(nearest_prot)

    def _conclude(self):
        self.results_distances = np.array(self.all_distances)
        self.results_nearest_resids = np.array(self.all_nearest_resids)


def get_lipid_protein_distances(u, prot_sel, lip_sel, step=1, use_com=False):
    """
    Wrapper function to execute continuous distance and nearest-residue calculation.
    """
    mode = "Center of Mass" if use_com else "Atomic"
    print(f"Calculating {mode} distances over {u.trajectory.n_frames} frames...")
    
    analysis = LipidProteinDistanceAnalysis(u, prot_sel, lip_sel, use_com=use_com)
    analysis.run(step=step, verbose=True)
    
    # We grab the resids we defined in __init__
    lipid_resids = analysis.lip_resids
    
    return lipid_resids, analysis.results_distances, analysis.results_nearest_resids