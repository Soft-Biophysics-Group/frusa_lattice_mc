"""
Vincent Ouazan-Reboul, 2025
Few functions to get the size of aggregates on any implemented lattice.
"""

import numpy as np
from pathlib import Path
import config as cfg
from geometry.load_geometry import geometry_from_model_file

from typing import TypeAlias
Aggregates: TypeAlias = list[set[int]]


def get_aggregates(
    struct_index: int = -1,
    struct_folder: str | Path = "",
    struct_file: str | Path = "",
    model_file: str | Path = cfg.default_mc_params_file,
) -> Aggregates:
    """
    Returns a list of sets. Each list member is a cluster, represented by the set of the
    sites containing its constituent particles.
    If struct_file is specified: overrides directly fetches results in
    struct_file.
    If struct_index is specified: load structure with index struct_index from folder
    struct_folder.
    If not, laods final structure from struct_folder.
    """
    site_orientations = cfg.load_structure(
        struct_index=struct_index, struct_folder=struct_folder, struct_file=struct_file
    )
    full_sites = cfg.get_full_sites(site_orientations)
    full_sites_set = set(full_sites)
    lattice, _ = geometry_from_model_file(model_file)

    visited_sites = set()
    first_shell = set()
    to_visit = set()
    all_clusters = []

    for site in full_sites:
        if site not in visited_sites:
            to_visit.add(site)
            cluster = set()
            # We build each cluster iteratively, shell by shell, by going through the neighbours
            # of each site and adding them to the cluster until we run out of them
            while len(to_visit) > 0:
                # Build first shell of the cluster
                site_to_visit = to_visit.pop()
                cluster.add(site_to_visit)
                visited_sites.add(site_to_visit)
                neighbour_sites = lattice.get_neighbour_sites(
                    site_to_visit
                )
                full_neighbour_sites = full_sites_set.intersection(set(neighbour_sites))
                # This is a set substraction: only keeps non-visited sites in the full neighbour
                # sites
                neighbours_to_visit = full_neighbour_sites - visited_sites
                # Updates the to_visit set with its intersection with neighbours_to_visit
                to_visit |= neighbours_to_visit
            all_clusters.append(cluster)

    return all_clusters

def get_aggregate_sizes(aggregates: Aggregates) -> list[int]:
    """
    Returns the sizes of all aggregates contained in Aggregates
    """
    return [len(agg) for agg in aggregates]
