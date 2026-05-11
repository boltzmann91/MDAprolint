import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import distance_array

class InteractingLipidsAnalysis(AnalysisBase):
    """
    Custom MDAnalysis tool to find interacting lipids and their nearest protein residue.
    """
    def __init__(self, universe, prot_sel, lip_sel, cutoff, **kwargs):
        super().__init__(universe.trajectory, **kwargs)
        self.prot_atoms = universe.select_atoms(prot_sel)
        self.lip_atoms = universe.select_atoms(lip_sel)
        self.cutoff = cutoff
        
        # Pre-extract resids for BOTH lipids and protein
        self.lip_resids = self.lip_atoms.resids
        self.prot_resids = self.prot_atoms.resids

    def _prepare(self):
        self.interactions_per_frame = []

    def _single_frame(self):
        # dist_matrix shape: (n_protein_atoms, n_lipid_atoms)
        dist_matrix = distance_array(
            self.prot_atoms.positions,
            self.lip_atoms.positions,
            box=self._ts.dimensions,
            backend='OpenMP'
        )
        
        # 1. Find the shortest distance for each lipid atom to ANY protein atom
        min_dists = dist_matrix.min(axis=0)
        
        # 2. Find the INDEX of that closest protein atom
        min_prot_indices = dist_matrix.argmin(axis=0)
        
        # 3. Create a True/False mask for lipid atoms within the cutoff
        valid_mask = min_dists <= self.cutoff
        
        # Extract the actual data only for atoms within the cutoff
        valid_lip_resids = self.lip_resids[valid_mask]
        valid_prot_resids = self.prot_resids[min_prot_indices[valid_mask]]
        valid_dists = min_dists[valid_mask]
        
        # 4. Map lipid_resid -> nearest_protein_resid
        interacting_dict = {}
        for l_res, p_res, dist in zip(valid_lip_resids, valid_prot_resids, valid_dists):
            if l_res not in interacting_dict:
                interacting_dict[l_res] = (p_res, dist)
            else:
                # If a lipid molecule has multiple atoms within the cutoff, 
                # strictly keep the protein residue that is absolutely closest
                if dist < interacting_dict[l_res][1]:
                    interacting_dict[l_res] = (p_res, dist)
                    
        # Strip the distance, keeping only the protein resid
        # Final result looks like: {105: 42, 106: 8, 202: 15}
        final_interactions = {k: v[0] for k, v in interacting_dict.items()}
        
        self.interactions_per_frame.append({
            'frame': self._ts.frame,
            'lipids': final_interactions
        })

    def _conclude(self):
        self.results = self.interactions_per_frame

def get_interacting_lipids(u, prot_sel, lip_sel, cutoff, step=1):
    """
    Finds interacting lipids and their primary binding site residue on the protein.
    """
    print(f"Finding interactions over {u.trajectory.n_frames} frames...")
    analysis = InteractingLipidsAnalysis(u, prot_sel, lip_sel, cutoff)
    analysis.run(step=step, verbose=True)
    return analysis.results
    