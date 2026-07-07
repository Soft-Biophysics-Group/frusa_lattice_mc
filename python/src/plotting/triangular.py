# pyright: basic
"""2D plotting of triangular-lattice (hexagonal) particles.

The lattice-independent machinery lives in `plotting_utils.ParticleRepresentation2D`; this module
only specifies the triangular-lattice details.
"""
import numpy as np
from .plotting_utils import ParticleRepresentation2D, ARROW_COLORS

# colormap of camembert contacts: forbidden contacts are in red, crystal in cyan, line in blue,
# nothing in orange
CAMEMBERT_CONTACTS_CMAP = ["red", "cyan", "blue", "orange"]
CAMEMBERT_CONTACTS = np.zeros((6, 6), dtype=int)
for i in range(6):
    CAMEMBERT_CONTACTS[i, (i + 3) % 6] = 1
CAMEMBERT_CONTACTS[0, 2] = 2
CAMEMBERT_CONTACTS[2, 0] = 2
CAMEMBERT_CONTACTS[1, 5] = 2
CAMEMBERT_CONTACTS[5, 1] = 2


class TriangularParticleRepresentation(ParticleRepresentation2D):
    """Plots hexagonal particles tiling a triangular lattice (6 faces, 6 neighbours)."""

    lattice_name = "triangular"
    n_faces = 6
    colors = [
        "bf9c76ff",
        "cf938dff",
        "b996c1ff",
        "7fa5d3ff",
        "67b0b0ff",
        "92ab7dff",
    ]


# Backwards-compatible alias
ParticleRepresentation = TriangularParticleRepresentation

__all__ = [
    "TriangularParticleRepresentation",
    "ParticleRepresentation",
    "ARROW_COLORS",
    "CAMEMBERT_CONTACTS",
    "CAMEMBERT_CONTACTS_CMAP",
]
