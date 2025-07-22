"""Vincent Ouazan-Reboul, 2025/07/22

Functions to get the lattice state associated to the figure shown in a .blend file
"""

from pathlib import Path
import bpy
from geometry import LatticeGeometry, ParticleGeometry
import numpy as np
from numpy.typing import NDArray
from scipy.spatial.transform import Rotation as R

file_path = "../../../../3dFigures/04_longerer_run_0.blend"

def get_objs_in_collection_locs_rots(
    filepath: Path | str, collection_name: str = "Particles"
):

    # Load .blend file
    with bpy.data.libraries.load(file_path, link = False) as (data_from, data_to):

        if collection_name not in data_from.collections:
            raise ValueError("Provided collection name not present in file")

        data_to.collections = [collection_name]

    coll = bpy.data.collections.get(collection_name)
    n_objs = len(coll.objects)

    all_locations = np.zeros((n_objs, 3))
    all_rotations_xyz_rad = np.zeros((n_objs, 3))

    for i, obj in enumerate(coll.objects):
        all_locations[i, :] = obj.location

        x_rot, y_rot, z_rot = obj.rotation_euler
        all_rotations_xyz_rad[i, :] = obj.rotation_euler

    return all_locations, all_rotations_xyz_rad

def auto_fit_lattice(
    particle_locations: NDArray[np.float_],
    lattice_name: str,
    make_dims_equal: bool = True,
) -> tuple[NDArray[np.float_], int, int, int]:
    # Make a ridiculously large lattice for our purposes
    infinite_lattice = LatticeGeometry.from_lattice_name(
        lattice_name, np.inf, np.inf, np.inf
    )
    # Convert all cartesian coordinates to infinite lattice to find span of lattice that fits
    # our structure snugly
    lattice_locations = np.zeros_like(particle_locations)
    for i, location in enumerate(particle_locations):
        # Note that we don't use cartesian_to_lattice here because it would apply PBCs.
        # At this point some of our lattice coordinates may stil be negative, and so this would
        # lead to infinities in our coordinates. Instead, we only multiply by the inverse basis
        # matrix.
        lattice_locations[i, :] = infinite_lattice.cartesian_to_lattice_basis @ location

    # print("Lattice_locations")
    # print(lattice_locations)
    # Offset particles so that they all have positive and integer lattice coordinates
    min_coords_lattice = np.min(lattice_locations, axis=0)
    # print("Minimum coordinates before offset")
    # print(min_coords_lattice)
    lattice_locations -= min_coords_lattice
    # print("Locations after offset")
    # print(lattice_locations)
    # print("Minimum coordinates after offset")
    # print(np.min(lattice_locations, 0))

    # Add 1 to avoid any PBC shenanigans
    lx, ly, lz = np.max(lattice_locations, 0).astype(int) + 1
    if make_dims_equal:
        max_dim = np.max([lx, ly, lz])
        lx, ly, lz = max_dim, max_dim, max_dim

    return lattice_locations, lx, ly, lz


def get_sites_and_orientations(
    particle_locations_cart: NDArray[np.float_],
    particle_rotations_xyz_rad: NDArray[np.float_],
    lattice_name: str,
    auto_fit_lattice_flag: bool = True,
    make_dims_equal:bool = True,
    lx: int | None = None,
    ly: int | None = None,
    lz: int | None = None,
):
    if auto_fit_lattice_flag:
        particle_locations_lattice, lx, ly, lz = auto_fit_lattice(
            particle_locations_cart, lattice_name, make_dims_equal = make_dims_equal
        )
    elif lx is None or ly is None or lz is None:
        raise ValueError("At least one lattice size missing and auto-fit disabled")

    lattice = LatticeGeometry.from_lattice_name(lattice_name, lx, ly, lz)
    particle = ParticleGeometry.from_lattice_name(lattice_name)

    if not auto_fit_lattice_flag:
        particle_locations_lattice = np.array(
            [
                lattice.cartesian_to_lattice(*location)
                for location in particle_locations_cart
            ],
            dtype=int,
        )

    # Identify the lattice sites
    all_sites = []
    for loc_latt in particle_locations_lattice:
        all_sites.append(lattice.lattice_coords_to_lattice_site(*loc_latt))
    all_sites = np.array(all_sites, dtype=int)

    # Identify the particle orientations
    all_orientations = []
    for orientation in particle_rotations_xyz_rad:
        this_orientation_rotation = R.from_euler("xyz", orientation)
        all_orientations.append(
            particle.identify_orientation(this_orientation_rotation)
        )
    all_orientations = np.array(all_orientations)

    return all_sites, all_orientations, [lx, ly, lz]

def write_sites_orientations_from_blend(
    blend_file: str | Path,
    output_file: str | Path,
    lattice_name: str,
    collection_name: str = "Particles",
    make_dims_equal: bool = True,
    lattice_dimensions: NDArray[np.int_]
    | list[int]
    | tuple[int, int, int]
    | None = None,
):
    particle_locations, particle_rotations_xyz_rad = get_objs_in_collection_locs_rots(
        blend_file, collection_name
    )
    auto_fit_lattice_flag: bool
    if lattice_dimensions is not None:
        lx, ly, lz = lattice_dimensions
        auto_fit_lattice_flag = False
    else:
        lx, ly, lz = None, None, None
        auto_fit_lattice_flag = True
    full_sites, orientations, [lx, ly, lz] = get_sites_and_orientations(
        particle_locations,
        particle_rotations_xyz_rad,
        lattice_name,
        auto_fit_lattice_flag=auto_fit_lattice_flag,
        make_dims_equal = make_dims_equal,
        lx=lx,
        ly=ly,
        lz=lz,
    )
    n_sites = lx * ly * lz

    # Create a table encoding the state of the lattice similarly to our simulation results
    sites_orientations = np.zeros((2, n_sites), dtype = int)
    sites_orientations[1, :] -= 1

    for (site, orientation) in zip(full_sites, orientations):
        sites_orientations[1, site] = orientation

    np.savetxt(output_file, sites_orientations)

    return sites_orientations

# Below: test operations

# lattice_name = "fcc"
#
# all_locations, all_rotations_xyz_rad = get_objs_in_collection_locs_rots(file_path)
# print("Locations")
# print(all_locations)
# locations_on_lattice, lx, ly, lz = auto_fit_lattice(all_locations, lattice_name)
#
# print(locations_on_lattice)
# print(f"{lx}, {ly}, {lz}")
#
# sites, orientations, lattice_dims = get_sites_and_orientations(
#     locations_on_lattice,
#     all_rotations_xyz_rad,
#     "fcc",
#     auto_fit_lattice_flag=False,
#     lx=lx,
#     ly=ly,
#     lz=lz,
# )
#
# output_file = "test_out.dat"
# sites_orientations = write_sites_orientations_from_blend(
#     file_path, output_file, "fcc",# lattice_dimensions=[lx, ly, lz]
# )
#
# print(f"{lx}, {ly}, {lz}")
#
# print(sites)
# print(orientations)
# print(sites_orientations)
# print(np.sum(sites_orientations[1, :] == -1))
# print(np.sum(sites_orientations[1, :] != -1))
