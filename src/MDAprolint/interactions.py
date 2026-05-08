import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import distance_array

class InteractingLipidsAnalysis(AnalysisBase):
    """
    Custom MDAnalysis tool to find interacting lipids using OpenMP distance arrays.
    """
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
    """
    Wrapper function to execute the OpenMP interaction calculation.
    """
    print(f"Finding interactions over {u.trajectory.n_frames} frames...")
    analysis = InteractingLipidsAnalysis(u, prot_sel, lip_sel, cutoff)
    analysis.run(step=step, verbose=True)
    return analysis.results
