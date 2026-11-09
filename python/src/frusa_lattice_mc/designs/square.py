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
        """Sides of the solid rectangle with the same gyration tensor as the cluster,
        largest first. A stripe is infinitely long; a cluster wrapping along several
        axes is nan.

        A solid n-site row of lattice sites has positional variance (n^2 - 1) / 12, so
        each side is sqrt(12 lambda + 1) for an eigenvalue lambda of the gyration
        tensor: exact for filled rectangles, and the area matches the size for
        compact clusters. Rejected: the rectangle with the same area and outer
        perimeter (calculation AAH, 26-07-22). Clusters without notches have the
        perimeter of their bounding box, so missing corner sites only lower the area,
        which that rectangle absorbs by stretching: an n x n square missing k sites
        gets width n - sqrt(k).
        """
        spanned = self.percolates_along
        if spanned.sum() > 1:
            return (np.nan, np.nan)
        if spanned.any():
            # Unbounded where it wraps, but the width still follows from the area
            return (np.inf, self.size / self.box_dims[spanned].item())

        coords = np.array(list(self.particle_coords_no_pbc.values()), dtype=float)
        # Flat directions (z in a 2D system) carry no extent
        coords = coords[:, self.box_dims > 1]
        eigenvalues = np.linalg.eigvalsh(np.cov(coords.T, bias=True))
        # Roundoff can push the zero eigenvalue of a straight line slightly negative
        sides = np.sqrt(12 * np.clip(eigenvalues, 0, None) + 1)
        return (sides[-1], sides[-2])

    @property
    def length(self) -> float:
        return self.effective_dims[0]

    @property
    def width(self) -> float:
        return self.effective_dims[1]


class RectangularAssemblies(Clusters):
    """The aggregates of a state, i.e. its connected components."""

    _cluster_class = RectangularAssembly
