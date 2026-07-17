""" Vincent Ouazan-Reboul, 2026
Storing the contact map designs for square lattice particles.
"""

from ..contact_utils import ContactMapWrapper
from ..geometry.particle_geometry import SquareParticle

CRYSTAL_CONTACTS = SquareParticle().get_all_canonical_contacts([(i, i+2) for i in range(2)])
SQUARE_VORTEX_CONTACTS = SquareParticle().get_all_canonical_contacts([(0, 1), (0, 3)])
