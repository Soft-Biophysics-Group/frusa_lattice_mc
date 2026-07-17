"""Vincent Ouazan-Reboul, 2025
Storing the contact map designs for triangular lattice particles (i.e. hexagons).
"""

from ..contact_utils import ContactMapWrapper
from ..geometry.particle_geometry import TriangularParticle

CRYSTAL_CONTACTS = TriangularParticle().get_all_canonical_contacts(
    [(i, i + 3) for i in range(3)]
)
VORTEX_CAMEMBERT_CONTACTS = TriangularParticle().get_all_canonical_contacts(
    [(0, 4), (1, 5)]
)
PATTERN_BULK_CONTACTS = (
    VORTEX_CAMEMBERT_CONTACTS
    | TriangularParticle().get_all_canonical_contacts([(2, 2), (3, 3)])
)

# -------- Geometry functions ---------


def calc_ec_ed_ratio_vortex(target_radius: float):
    return (2 * target_radius**2 - 2 * target_radius - 1) / (3 * target_radius**2)


def calc_ec_ed_ratio_patterned_bulk(target_radius: float):
    return (2 * target_radius**2 - 2 * target_radius - 1) / (2 * target_radius**2)
