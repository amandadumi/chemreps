'''
Initializing of the main representation functions
'''
from . import utils
from . import bagger
from . import coulomb_matrix
from . import bag_of_bonds
from . import bat
from . import just_bonds
from . import dataset
try:
    import rdkit
    from . import fingerprints
except ImportError:
   print('Warning: Optional dependency rdkit not installed, the fingerprints representation will not be available')

