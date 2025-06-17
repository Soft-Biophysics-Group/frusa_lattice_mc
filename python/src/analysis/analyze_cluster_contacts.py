from geometry import LatticeGeometry, ParticleGeometry
from config import load_model_file, load_structure
import config as cfg
from pathlib import Path
from analysis.analyze_3d_size import get_aggregates, get_agg_sizes
from analysis.get_face_contacts import get_particle_face_face_contacts

contact_type = dict[frozenset[int], frozenset[int]]


def get_one_aggregate_contacts(
    aggregate: set[int],
    face_face_contacts: contact_type,
    lg: LatticeGeometry,
    pg: ParticleGeometry,
    contacts_to_check: list[frozenset[int]],
):

    checked_pairs:list[frozenset[int]] = []

    n_contacts = 0

    for site_1 in aggregate:
        neighbours_in_agg = set(lg.get_neighbour_sites(site_1)) & aggregate
        for site_2 in neighbours_in_agg:
            pair_set = frozenset((site_1, site_2))
            if pair_set not in checked_pairs:
                face_face_contact = face_face_contacts[pair_set]
                if face_face_contact in contacts_to_check:
                    n_contacts += 1

                checked_pairs.append(pair_set)

    return n_contacts

def get_all_aggregate_contacts(
    all_contact_types_to_check: list[list[tuple[int]]] | list[list[frozenset[int]]],
    struct_index: int = -1,
    struct_folder: str | Path = "",
    struct_file: str | Path = "",
    model_file: str | Path = cfg.default_mc_params_file,
):
    """Generates aggregates from simulation results and counts the number of certain specific
    contacts, grouped into classes.
    """

    # Get the simulation data
    sites_orientations = load_structure(struct_index, struct_folder, struct_file)
    full_sites = cfg.get_full_sites(sites_orientations)
    face_face_contacts = get_particle_face_face_contacts(
        struct_index, struct_folder, struct_file, model_file = model_file
    )
    lattice_geometry = LatticeGeometry.from_model_file(model_file)
    particle_geometry = ParticleGeometry.from_model_file(model_file)
    all_aggregates = get_aggregates(
        struct_index, struct_folder, struct_file, model_file=model_file
    )

    # Aggregate variable for results
    results: dict[str, list[int] | list[list[int]] | list[set[int]]] = {
        "sites": all_aggregates
    }

    # Check the number of contacts of each type for each cluster
    all_agg_contact_types = []
    for aggregate in all_aggregates:
        this_agg_contact_types = []
        for contact_type in all_contact_types_to_check:
            contacts_to_check_set = [
                particle_geometry.get_canonical_contact(*contact)
                for contact in contact_type
            ]
            n_contacts = get_one_aggregate_contacts(
                aggregate,
                face_face_contacts,
                lattice_geometry,
                particle_geometry,
                contacts_to_check_set,
            )
            this_agg_contact_types.append(n_contacts)
        all_agg_contact_types.append(this_agg_contact_types)
    results["n_contact_agg_type"] = all_agg_contact_types

    all_sizes = [len(agg) for agg in all_aggregates]
    results["sizes"] = all_sizes

    return results
