"""Vincent Ouazan-Reboul, 2025
Set of functions to find the configurations (using periodic boundary conditions) of a simulation
output which minimize the amount of periodic boundary conditions.
"""

import numpy as np
from numpy.typing import NDArray
from geometry import LatticeGeometry
import config as cfg
from tqdm import tqdm

def is_pbc_bond(bond:list[int] | NDArray[np.int_]) -> bool:
    return np.max(np.abs(bond)) > 1

def count_n_pbc_contacts(full_sites: NDArray[np.int_], lattice_geometry:LatticeGeometry) -> int:
    checked_pairs: list[set[int]] = []
    full_sites_set: set[int] = set(full_sites)

    n_pbc_bonds = 0

    for site_1 in full_sites_set:
        coords_of_1 = lattice_geometry.lattice_site_to_lattice_coords(site_1)
        neighbours_of_1 = set(lattice_geometry.get_neighbour_sites(site_1))
        full_neighbours_of_1 = (neighbours_of_1 & full_sites_set)
        for site_2 in full_neighbours_of_1:
            if {site_1, site_2} not in checked_pairs:
                checked_pairs.append({site_1, site_2})
                coords_of_2 = lattice_geometry.lattice_site_to_lattice_coords(site_2)
                naive_bond = coords_of_2 - coords_of_1
                n_pbc_bonds += is_pbc_bond(naive_bond)

    return n_pbc_bonds


def count_all_possible_n_pbc_contacts(
    site_orientations: NDArray[np.int_], lattice_geometry: LatticeGeometry
) -> NDArray[np.int_]:
    all_n_pbc = []
    n_sites = lattice_geometry.lx * lattice_geometry.ly * lattice_geometry.lz
    for site in tqdm(range(n_sites)):
        site_vector = lattice_geometry.lattice_site_to_lattice_coords(site)
        translated_site_orientations = lattice_geometry.apply_translation_to_config(
            site_orientations, site_vector
        )
        translated_full_sites = cfg.get_full_sites(translated_site_orientations)
        all_n_pbc.append(count_n_pbc_contacts(translated_full_sites, lattice_geometry))

    return np.array(all_n_pbc)


def find_best_translation(struct_file, model_file):
    site_orientations = cfg.load_structure(struct_file=struct_file)
    lattice_geometry = LatticeGeometry.from_model_file(model_file)
    all_pbc_counts = count_all_possible_n_pbc_contacts(
        site_orientations, lattice_geometry
    )

    lowest_pbc_number = np.min(all_pbc_counts)
    best_sites = np.where(all_pbc_counts == lowest_pbc_number)[0]

    best_translations = [
        lattice_geometry.lattice_site_to_lattice_coords(site) for site in best_sites
    ]

    return best_translations, lowest_pbc_number
