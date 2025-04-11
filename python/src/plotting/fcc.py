"""
Vincent Ouazan-Reboul, 2025
Functions to plot the results of simulations of FCC particles.
So far only supports one type of particles.
"""

from pathlib import Path

path_to_config = Path(__file__).parent.parent

ALL_CONTACTS = [(i, j) for i in range(24) for j in range(24)]
# DEFAULT_MATERIAL = (14, 0, 255, 0.1)
DEFAULT_MATERIAL = (14, 0, 255, 1.0)

src_plotting_dir = Path(__file__).parent

path_to_numbered_rhombic = (
    src_plotting_dir / "assets/oneRhombicDodecahedron/one_rhombic_numbered.obj"
)
path_to_one_axis_rhombic = (
    src_plotting_dir / "assets/oneRhombicDodecahedron/one_rhombic_one_axis.obj"
)

def plot_boundary(*args, **kwargs):
    return
