import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import distance_array

class InteractingLipidsAnalysis(AnalysisBase):
    """
    Custom MDAnalysis tool to find interacting lipids using Dual-Cutoff logic.
    Supports atomic positions or Centers of Mass (COM).
    """
    def __init__(self, universe, prot_sel, lip_sel, lower_cutoff, upper_cutoff, use_com=False, **kwargs):
        super().__init__(universe.trajectory, **kwargs)
        self.prot_atoms = universe.select_atoms(prot_sel)
        self.lip_atoms = universe.select_atoms(lip_sel)
        
        self.lower_cutoff = lower_cutoff
        self.upper_cutoff = upper_cutoff
        self.use_com = use_com
        
        # Pre-extract resids for mapping
        if self.use_com:
            self.lip_resids = self.lip_atoms.residues.resids
        else:
            self.lip_resids = self.lip_atoms.resids
            
        self.prot_resids = self.prot_atoms.resids

    def _prepare(self):
        self.interactions_per_frame = []
        self.active_bound_lipids = {} 

    def _single_frame(self):
        # 1. Get Lipid Positions
        if self.use_com:
            lip_positions = self.lip_atoms.center_of_mass(compound='residues')
        else:
            lip_positions = self.lip_atoms.positions

        # 2. Distance Matrix
        dist_matrix = distance_array(
            self.prot_atoms.positions,
            lip_positions,
            box=self._ts.dimensions,
            backend='OpenMP'
        )
        
        min_dists = dist_matrix.min(axis=0)
        min_prot_indices = dist_matrix.argmin(axis=0)
        
        tight_mask = min_dists <= self.lower_cutoff
        loose_mask = min_dists <= self.upper_cutoff
        
        # 3. Map lipid_resid -> nearest protein_resid for ALL loose atoms/COMs
        loose_dict = {}
        for l_res, p_res, dist in zip(self.lip_resids[loose_mask], 
                                      self.prot_resids[min_prot_indices[loose_mask]], 
                                      min_dists[loose_mask]):
            if l_res not in loose_dict or dist < loose_dict[l_res][1]:
                loose_dict[l_res] = (p_res, dist)
                
        tight_lipid_resids = set(self.lip_resids[tight_mask])
        
        # --- DUAL CUTOFF STATE MACHINE ---
        for l_res in tight_lipid_resids:
            self.active_bound_lipids[l_res] = loose_dict[l_res][0]
            
        active_keys = list(self.active_bound_lipids.keys())
        for l_res in active_keys:
            if l_res in loose_dict:
                self.active_bound_lipids[l_res] = loose_dict[l_res][0]
            else:
                del self.active_bound_lipids[l_res]
                
        self.interactions_per_frame.append({
            'frame': self._ts.frame,
            'lipids': dict(self.active_bound_lipids) 
        })

    def _conclude(self):
        self.results = self.interactions_per_frame


def get_interacting_lipids(u, prot_sel, lip_sel, lower_cutoff, upper_cutoff, step=1, use_com=False):
    """
    Finds interacting lipids using a mathematically rigorous dual-cutoff approach.
    """
    mode = "Center of Mass" if use_com else "Atomic"
    print(f"Finding interactions over {u.trajectory.n_frames} frames ({mode} | Dual Cutoffs: {lower_cutoff:.2f} / {upper_cutoff:.2f} Å)...")
    
    analysis = InteractingLipidsAnalysis(u, prot_sel, lip_sel, lower_cutoff, upper_cutoff, use_com=use_com)
    analysis.run(step=step, verbose=True)
    return analysis.results