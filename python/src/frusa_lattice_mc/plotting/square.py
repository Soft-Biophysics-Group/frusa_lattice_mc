# pyright: basic
"""2D plotting of square-lattice particles.

The lattice-independent machinery lives in `plotting_utils.ParticleRepresentation2D`; this module
only specifies the square-lattice details.
"""
from .plotting_utils import ParticleRepresentation2D, ARROW_COLORS


class SquareParticleRepresentation(ParticleRepresentation2D):
    """Plots square particles tiling a square lattice (4 faces, 4 neighbours)."""

    lattice_name = "square"
    n_faces = 4
    colors = [
        "bf9c76ff",
        "b996c1ff",
        "7fa5d3ff",
        "92ab7dff",
    ]


# Alias mirroring plotting.triangular for a uniform per-lattice interface
ParticleRepresentation = SquareParticleRepresentation

__all__ = ["SquareParticleRepresentation", "ParticleRepresentation", "ARROW_COLORS"]
