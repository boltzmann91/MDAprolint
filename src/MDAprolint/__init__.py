# Import the specific functions you want to expose to the user
from .core import get_interacting_lipids
from .core import calculate_residence_times
from .core import get_lipid_protein_distances

# Optional: define __all__ to restrict what is imported with 'from lipid_interactions import *'
__all__ = [
    "get_interacting_lipids",
    "calculate_residence_times",
    "get_lipid_protein_distances"
]
