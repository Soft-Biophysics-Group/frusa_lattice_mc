""" Vincent Ouazan-Rebo, 2025
Storing the contact map designs for triangular lattice particles (i.e. hexagons).
"""

from geometry.cubic import CubicParticle, CubicGeometry
import config as cfg
from contact_utils import ContactMapWrapper

CRYSTAL_CONTACTS = [(i, i+3%6) for i in range(6)]

def set_crystal_contacts(crystal_e:float, cmap: ContactMapWrapper):
    for contact in CRYSTAL_CONTACTS:
        cmap[*contact] = crystal_e
    return

def set_vortex_camembert_contacts(defect_e, cmap:ContactMapWrapper):

