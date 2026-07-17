"""
Vincent Ouazan-Reboul, 2025
Few functions to get the size of aggregates on any implemented lattice.
"""

from pathlib import Path
from .. import config as cfg
from ..geometry import LatticeGeometry
from ..lattice_state import LatticeState

from typing import TypeAlias

Aggregates: TypeAlias = list[set[int]]

class Cluster:
    particles: list[int]
    size: int
    face_face_contacts: dict[frozenset[int], frozenset[int]]


def get_aggregates(
    lattice_state: LatticeState,
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

    visited_sites = set()
    to_visit = set()
    all_clusters = []

    for site in lattice_state.full_sites:
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
                neighbour_sites = lattice_state.lattice.get_neighbour_sites(site_to_visit)
                full_neighbour_sites = set(lattice_state.full_sites).intersection(
                    set(neighbour_sites)
                )
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


def get_crystalline_domains(
    struct_index: int | None = None,
    struct_folder: str | Path | None = None,
    struct_file: str | Path | None = None,
    model_file: str | Path = cfg.default_mc_params_file,
) -> Aggregates:
    """
    Returns a list of sets. Each member is a crystalline domain, i.e. a set of connected particles with the same orientation.
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
    lattice = LatticeGeometry.from_model_file(model_file)

    visited_sites = set()
    to_visit = set()
    all_crystalline_components = []

    for site in full_sites:
        if site not in visited_sites:
            to_visit.add(site)
            cryst_comp = set()
            ref_orientation = site_orientations[1, site]
            # We build each cluster iteratively, shell by shell, by going through the neighbours
            # of each site and adding them to the cluster until we run out of them
            while len(to_visit) > 0:
                # Build first shell of the cluster
                site_to_visit = to_visit.pop()
                orientation_site_to_visit = site_orientations[1, site_to_visit]
                if orientation_site_to_visit == ref_orientation:
                    cryst_comp.add(site_to_visit)
                    visited_sites.add(site_to_visit)
                    neighbour_sites = lattice.get_neighbour_sites(site_to_visit)
                    full_neighbour_sites = full_sites_set.intersection(
                        set(neighbour_sites)
                    )
                    # This is a set substraction: only keeps non-visited sites in the full neighbour
                    # sites
                    neighbours_to_visit = full_neighbour_sites - visited_sites
                    # Updates the to_visit set with its intersection with neighbours_to_visit
                    to_visit |= neighbours_to_visit
            all_crystalline_components.append(cryst_comp)

    return all_crystalline_components


def get_n_interfaces_w_ext(
    aggregates: Aggregates,
    model_file: str | Path = cfg.default_mc_params_file,
) -> list[int]:
    """Returns the number of interfaces between each aggregate and the outside."""
    n_interfaces = []
    lattice = LatticeGeometry.from_model_file(model_file)
    all_full_sites = set()

    for agg in aggregates:
        all_full_sites |= agg

    for agg in aggregates:
        this_n_interfaces = 0
        for site in agg:
            neighbours = lattice.get_neighbour_sites(site)
            for neighbour in neighbours:
                if neighbour not in all_full_sites:
                    this_n_interfaces += 1
        n_interfaces.append(this_n_interfaces)

    return n_interfaces
