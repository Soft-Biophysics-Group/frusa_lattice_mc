"""Vincent Ouazan-Reboul, 2025
Storing the contact map designs for triangular lattice particles (i.e. hexagons).
"""

from ..geometry.particle_geometry import TriangularParticle
import numpy as np
from functools import cached_property
from ..analysis.clusters import Cluster, Clusters

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

def calc_radius_of_ec_ed_ratio_vortex(ec_ed:float):
    return (1 + np.sqrt(3 * (1 - ec_ed))) / (2 - 3 * ec_ed)

def calc_vortex_radius_of_size(n_particles: int) -> float:
    return (-1 + np.sqrt(1 + 4 * n_particles / 3)) / 2

def calc_vortex_size_of_radius(radius: float) -> float:
    return 3 * radius * (radius + 1)


def calc_ec_ed_ratio_patterned_bulk(target_radius: float):
    return (2 * target_radius**2 - 2 * target_radius - 1) / (2 * target_radius**2)


# ---------- Cluster specialization to our aggregates of interest ----------
class VortexAssembly(Cluster):
    @cached_property
    def radius(self) -> float:
        perimeter = self.outer_perimeter
        area = self.size
        # Monomer is radius 0
        if area == 1 and perimeter == 6:
            return 0
        # If discriminant is negative, return 0 (I don't think this ever happens)
        discriminant = 9 - 12 * (perimeter / 2 - area - 3)
        if discriminant < 0:
            return 0
        else:
            return (3 + np.sqrt(discriminant)) / 6

    @cached_property
    def total_branch_length(self) -> float:
        if self.radius == 0:
            return 0
        return self.outer_perimeter / 2 - 3 * (2 * self.radius + 1)


class VortexAssemblies(Clusters):
    _cluster_class = VortexAssembly

    def __init__(self, state: LatticeState):
        super().__init__(state, find_connected_site_sets(state, lambda x, y: True))

class PatternedBulkTriangle(Cluster):
    @cached_property
    def sector_size(self) -> float:
        return 6 * self.size


class PatternedBulks(Clusters):
    _cluster_class = PatternedBulkTriangle

    def __init__(self, state: LatticeState):
        orientations = state.orientations
        super().__init__(
            state,
            find_connected_site_sets(
                state, lambda x, y: orientations[x] == orientations[y]
            ),
        )
