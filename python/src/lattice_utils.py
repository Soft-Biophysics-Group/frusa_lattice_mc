# Vincent Ouazan-Reboul
# 2024/12/30
# Python utilities we will be using to handle all lattices in our data analysis

import numpy as np

# def lattice_site_to_lattice_coords_2d(site_index: int, lx: int):
#     y_lattice = site_index // lx
#     x_lattice = site_index - lx * y_lattice
#     return np.array([x_lattice, y_lattice])
#
#
# def lattice_coords_to_lattice_site_2d(x_lattice: int, y_lattice: int, lx: int):
#     return x_lattice + y_lattice * lx
#
def get_full_sites_characteristics(results):
    """
    Taking in a results array which has the same format as the one returned by the c++ program,
    returns a 3D array with the occupied site as the first column, the particle type as the
    second, and the particle orientation as the third.
    """
    full_sites = []
    for site, (ptype, orientation) in enumerate(results.T):
        if orientation != -1:
            full_sites.append([site, ptype, orientation])
    return np.vstack(full_sites)
