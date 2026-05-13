import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import distance_array

class InteractingLipidsAnalysis(AnalysisBase):
    """
    Custom MDAnalysis tool to find interacting lipids using Dual-Cutoff logic.
    """
    def __init__(self, universe, prot_sel, lip_sel, lower_cutoff, upper_cutoff, **kwargs):
        super().__init__(universe.trajectory, **kwargs)
        self.prot_atoms = universe.select_atoms(prot_sel)
        self.lip_atoms = universe.select_atoms(lip_sel)
        
        self.lower_cutoff = lower_cutoff
        self.upper_cutoff = upper_cutoff
        
        self.lip_resids = self.lip_atoms.resids
        self.prot_resids = self.prot_atoms.resids

    def _prepare(self):
        self.interactions_per_frame = []
        
        # --- NEW: Dual Cutoff Memory ---
        # Tracks lipids currently bound: {lipid_resid: nearest_prot_resid}
        self.active_bound_lipids = {} 

    def _single_frame(self):
        dist_matrix = distance_array(
            self.prot_atoms.positions,
            self.lip_atoms.positions,
            box=self._ts.dimensions,
            backend='OpenMP'
        )
        
        min_dists = dist_matrix.min(axis=0)
        min_prot_indices = dist_matrix.argmin(axis=0)
        
        # 1. Identify atoms within cutoffs
        tight_mask = min_dists <= self.lower_cutoff
        loose_mask = min_dists <= self.upper_cutoff
        
        # 2. Map lipid_resid -> nearest protein_resid for ALL loose atoms
        loose_dict = {}
        for l_res, p_res, dist in zip(self.lip_resids[loose_mask], 
                                      self.prot_resids[min_prot_indices[loose_mask]], 
                                      min_dists[loose_mask]):
            if l_res not in loose_dict or dist < loose_dict[l_res][1]:
                loose_dict[l_res] = (p_res, dist)
                
        # 3. Get unique lipid resids that crossed the TIGHT (lower) cutoff
        tight_lipid_resids = set(self.lip_resids[tight_mask])
        
        # --- DUAL CUTOFF STATE MACHINE ---
        
        # Rule A: Any lipid that crosses the lower cutoff becomes ACTIVE
        for l_res in tight_lipid_resids:
            # Add to active memory and record its nearest protein residue
            self.active_bound_lipids[l_res] = loose_dict[l_res][0]
            
        # Rule B: Keep ACTIVE lipids ONLY if they are still within the upper cutoff
        # (We use list() to safely iterate while deleting dictionary keys)
        active_keys = list(self.active_bound_lipids.keys())
        for l_res in active_keys:
            if l_res in loose_dict:
                # Update the binding site (lipids can slide around the pocket while bound!)
                self.active_bound_lipids[l_res] = loose_dict[l_res][0]
            else:
                # The lipid has drifted past the upper cutoff. Evict it!
                del self.active_bound_lipids[l_res]
                
        # Save a COPY of the currently active lipids for this frame
        self.interactions_per_frame.append({
            'frame': self._ts.frame,
            'lipids': dict(self.active_bound_lipids) 
        })

    def _conclude(self):
        self.results = self.interactions_per_frame


def get_interacting_lipids(u, prot_sel, lip_sel, lower_cutoff, upper_cutoff, step=1):
    """
    Finds interacting lipids using a mathematically rigorous dual-cutoff approach.
    """
    print(f"Finding interactions over {u.trajectory.n_frames} frames (Dual Cutoffs: {lower_cutoff:.2f} / {upper_cutoff:.2f} Å)...")
    analysis = InteractingLipidsAnalysis(u, prot_sel, lip_sel, lower_cutoff, upper_cutoff)
    analysis.run(step=step, verbose=True)
    return analysis.results