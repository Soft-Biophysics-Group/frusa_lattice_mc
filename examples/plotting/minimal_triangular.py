"""Minimal example: plot a triangular-lattice (hexagonal particle) configuration.

Run with:  uv run python examples/plotting/minimal_triangular.py
"""
from pathlib import Path
import numpy as np
from frusa_lattice_mc.plotting.triangular import TriangularParticleRepresentation, CAMEMBERT_CONTACTS

LX, LY = 8, 8
OUT = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)

rep = TriangularParticleRepresentation(lx=LX, ly=LY)

# A structure is a (2, lx*ly) int array: row 0 = particle type, row 1 = orientation
# (-1 marks an empty site). Here we fill a small block with varied orientations.
structure = np.full((2, LX * LY), -1, dtype=int)
for x in range(2, 6):
    for y in range(2, 6):
        site = rep.lattice.lattice_coords_to_lattice_site(x, y)
        structure[:, site] = (0, (x + y) % rep.n_faces)

struct_file = OUT / "triangular_structure.dat"
np.savetxt(struct_file, structure, fmt="%d")

# Outlines + orientation arrows, then a contact overlay colouring the camembert contacts.
fig, ax = rep.plot_results_arrows(results_file=struct_file)
contact_colors = np.array(["black", "cyan", "blue", "orange"])[CAMEMBERT_CONTACTS]
rep.plot_contacts(ax, contact_colors, results_file=struct_file)
ax.set_aspect("equal")

fig.savefig(OUT / "minimal_triangular.png", dpi=150)
print("saved", OUT / "minimal_triangular.png")
