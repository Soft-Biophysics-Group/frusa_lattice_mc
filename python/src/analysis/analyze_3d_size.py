"""
Vincent Ouazan-Reboul, 2025
Few functions to get the size of aggregates on 3D lattices.
"""

import numpy as np
from pathlib import Path
import config as cfg
from geometry import LatticeGeometry


def get_aggregates(
    struct_index: int = -1,
    struct_folder: str | Path | None = None,
    struct_file: str | Path | None = None,
    model_file: str | Path = cfg.default_mc_params_file,
) -> list[set[int]]:
    """
    Returns a list of sets. Each list member is a cluster, represented by the set of the
    sites containing its constitueent particles.
    If struct_file is specified: overrides directly fetches results in
    struct_file.
    If struct_index is specified: load structure with index struct_index from folder
    struct_folder.
    If not, laods final structure from struct_folder.
    Returned structure is a dimensional numpy array, with line 1 corresponding to particle type
    and line 2 to particle orientation.
    orientation -1 always corresponds to an empty site, irrespective of particle type.
    """
    site_orientations = cfg.load_structure(
        struct_index=struct_index, struct_folder=struct_folder, struct_file=struct_file
    )
    full_sites = cfg.get_full_sites(site_orientations)
    full_sites_set = set(full_sites)
    lattice = LatticeGeometry.from_model_file(model_file)

    visited_sites = set()
    first_shell = set()
    to_visit = set()
    all_clusters = []

    # Build clusters neighbour by neighbour
    for site in full_sites:
        # Unvisited site starts a new cluster
        if site not in visited_sites:
            to_visit.add(site)
            cluster = set()
            # We have built the whole cluster if we run out of sites to visit
            while len(to_visit) > 0:
                site_to_visit = to_visit.pop()
                # Using a set guarantees we don't have more than one copy of each site
                cluster.add(site_to_visit)
                visited_sites.add(site_to_visit)
                neighbour_sites = lattice.get_neighbour_sites(site_to_visit)
                # Select only neighbouring sites that are full
                full_neighbour_sites = full_sites_set.intersection(set(neighbour_sites))
                # And visit only the unvisited ones
                neighbours_to_visit = full_neighbour_sites - visited_sites
                # Shorthand for set = union(set, other_set)
                to_visit |= neighbours_to_visit
            all_clusters.append(cluster)

    return all_clusters

def get_agg_sizes(aggregates: list[set[int]]) -> dict[int, int]:
    all_sizes = {}
    for agg in aggregates:
        agg_size = len(agg)
        if agg_size in all_sizes.keys():
            all_sizes[agg_size] += 1
        else:
            all_sizes[agg_size] = 1
    return all_sizes
