""" Vincent Ouazan-Reboul, 2026
Storing the contact map designs for square lattice particles.
"""

import numpy as np
from functools import cached_property

from ..contact_utils import ContactMapWrapper
from ..geometry.particle_geometry import SquareParticle
from ..analysis.clusters import Cluster, Clusters

CRYSTAL_CONTACTS = SquareParticle().get_all_canonical_contacts([(i, i+2) for i in range(2)])
SQUARE_VORTEX_CONTACTS = SquareParticle().get_all_canonical_contacts([(0, 1), (0, 3)])


# ---------- Cluster specialization to our aggregates of interest ----------
class RectangularAssembly(Cluster):
    @cached_property
    def effective_dims(self) -> tuple[float, float]:
        """Dimensions of the rectangle with the same area and perimeter as the cluster,
        largest first. See calculation AAH, 26-07-22 for derivation. A stripe is
        infinitely long, a shape the rectangle cannot fit is nan."""
        spanned = self.percolates_along
        if spanned.sum() > 1:
            return (np.nan, np.nan)
        if spanned.any():
            # Unbounded where it wraps, but the width still follows from the area
            return (np.inf, self.size / self.box_dims[spanned].item())

        area, perimeter = self.size, self.outer_perimeter
        # Ragged clusters can carry too much perimeter for their area
        discriminant = 1 - 16 * area / perimeter**2
        if discriminant < 0:
            return (np.nan, np.nan)
        root = np.sqrt(discriminant)
        return (perimeter / 4 * (1 + root), perimeter / 4 * (1 - root))

    @property
    def length(self) -> float:
        return self.effective_dims[0]

    @property
    def width(self) -> float:
        return self.effective_dims[1]


class RectangularAssemblies(Clusters):
    """The aggregates of a state, i.e. its connected components."""

    _cluster_class = RectangularAssembly
