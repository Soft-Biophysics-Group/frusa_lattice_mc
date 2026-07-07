"""Minimal example: plot a square-lattice configuration, selecting the lattice
automatically from a model file.

Run with:  uv run python examples/plotting/minimal_square.py
"""
import json
from pathlib import Path
import numpy as np
from plotting.plotting_utils import ParticleRepresentation2D

LX, LY = 8, 8
OUT = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)

# Only the model file names the lattice; the plotting code picks the matching
# representation (here SquareParticleRepresentation) on its own.
model_file = OUT / "square_model.json"
json.dump({"lattice_name": "square", "lx": LX, "ly": LY, "lz": 1}, open(model_file, "w"))
rep = ParticleRepresentation2D.from_model_file(model_file)

# A structure is a (2, lx*ly) int array: row 0 = particle type, row 1 = orientation
# (-1 marks an empty site). Here we fill a small block with varied orientations.
structure = np.full((2, LX * LY), -1, dtype=int)
for x in range(2, 6):
    for y in range(2, 6):
        site = rep.lattice.lattice_coords_to_lattice_site(x, y)
        structure[:, site] = (0, (x + y) % rep.n_faces)

struct_file = OUT / "square_structure.dat"
np.savetxt(struct_file, structure, fmt="%d")

fig, ax = rep.plot_results_arrows(results_file=struct_file)
ax.set_aspect("equal")

fig.savefig(OUT / "minimal_square.png", dpi=150)
print("saved", OUT / "minimal_square.png")
