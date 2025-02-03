"""
Vincent Ouazan-Reboul, 2025
Few functions to get the size of aggregates on 3D lattices.
"""

import numpy as np
from geometry.geometry import lattice_site_to_lattice_coords_3d, lattice_coords_to_lattice_site_3d
from pathlib import Path
import config as cfg

def get_neighbouring_particles(site, full_sites, lx, ly, lz):
    """
    Returns a set of the sites neighbouring `site` and containing a particle
    - `site` is the address of a site on the 3d lattice
    - `full_sites` is a list of full sites, typically obtained through `config.get_full_sites`
      function
    - lx, ly, lz are the lattice dimensions, which you can obtain using the
      `config.load_model_file` function
    """
    neighbours = set()

    x, y, z = lattice_site_to_lattice_coords_3d(site, lx, ly)
    xp1 = ((x + 1)%lx, y, z)
    yp1 = (x, (y + 1)%ly, z)
    zp1 = (x, y, (z + 1)%lz)
    xm1 = (x-1 if x > 0 else lx-1, y, z)
    ym1 = (x, y-1 if y > 0 else ly-1, z)
    zm1 = (x, y, z-1 if z > 0 else lz-1)

    for neigh_coords in [xp1, yp1, zp1, xm1, ym1, zm1]:
        xn, yn, zn = neigh_coords
        neigh_site = lattice_coords_to_lattice_site_3d(xn, yn, zn, lx, ly)
        if neigh_site in full_sites:
            neighbours.add(neigh_site)

    return neighbours


def get_aggregates(
    struct_index: int = -1,
    struct_folder: str | Path = "",
    struct_file: str | Path = "",
    model_file: str | Path = cfg.default_mc_params_file,
):
    """
    Returns a list of sets. Each list member is a cluster, represented by the set of the
    sites containing its constitueent particles.
    """
    site_orientations = cfg.load_structure(
        struct_index=struct_index, struct_folder=struct_folder, struct_file=struct_file
    )
    model_params = cfg.load_model_file(model_file)
    full_sites = cfg.get_full_sites(site_orientations)
    lx = model_params["lx"]
    ly = model_params["ly"]
    lz = model_params["lz"]

    visited_sites = set()
    first_shell = set()
    to_visit = set()
    all_clusters = []

    for site in full_sites:
        if site not in visited_sites:
            to_visit.add(site)
            cluster = set()
            while len(to_visit) > 0:
                # Build first shell of the cluster
                site_to_visit = to_visit.pop()
                cluster.add(site_to_visit)
                visited_sites.add(site_to_visit)
                neighbour_sites = get_neighbouring_particles(site_to_visit, full_sites, lx, ly, lz)
                neighbours_to_visit = neighbour_sites - visited_sites
                to_visit |= neighbours_to_visit
            all_clusters.append(cluster)

    return all_clusters
