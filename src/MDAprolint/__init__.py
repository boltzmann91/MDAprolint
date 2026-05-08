# Import the functions from their specific modules
from .interactions import get_interacting_lipids_openmp
from .kinetics import calculate_residence_times
from .distances import get_lipid_protein_distances

# Expose them to the user
__all__ = [
    "get_interacting_lipids_openmp",
    "calculate_residence_times",
    "get_lipid_protein_distances"
]

