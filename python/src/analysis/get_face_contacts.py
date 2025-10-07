from geometry.lattice_geometry import LatticeGeometry
from geometry.particle_geometry import ParticleGeometry
from pathlib import Path
import config as cfg


def get_particle_face_face_contacts(
    struct_index: int = -1,
    struct_folder: str | Path = "",
    struct_file: str | Path = "",
    model_file: str | Path = cfg.default_mc_params_file,
):
    """Returns all the particle face-face contacts as a dict for data analysis.

    Args:
    If struct_file is specified: overrides directly fetches results in
    struct_file.
    If struct_index is specified: load structure with index struct_index from folder
    struct_folder.
    If not, laods final structure from struct_folder.
    Returned structure is a dimensional numpy array, with line 1 corresponding to particle type
    and line 2 to particle orientation.

    Returns:
    face_face_contacts, a dict mapping sets of two full sites to a list of sets. Each set
    corresponds to one set of equivalent faces in contact of the two occupying particles
    (obtained by rotating the two particles around the face in which they contact). Keys are
    only the sets corresponding to two full sites!
    """

    site_orientations = cfg.load_structure(
        struct_index=struct_index, struct_folder=struct_folder, struct_file=struct_file
    )
    full_sites = cfg.get_full_sites(site_orientations)
    orientations = site_orientations[1, :]
    lattice = LatticeGeometry.from_model_file(model_file)
    particle = ParticleGeometry.from_model_file(model_file)

    all_contacts = {}

    for i, site_1 in enumerate(full_sites):
        orientation_1 = orientations[site_1]
        neighbours_of_1 = lattice.get_neighbour_sites(site_1)
        for neighbour in neighbours_of_1:
            if neighbour in full_sites:
                site_2 = neighbour
                particles_set = frozenset((site_1, site_2))
                if particles_set not in all_contacts.keys():
                    orientation_2 = orientations[site_2]
                    face_1, face_2, bond = lattice.get_faces_in_contact_and_bond(
                        site_1, orientation_1, site_2, orientation_2
                    )
                    canonical_contact = particle.get_canonical_contact(face_1, face_2)
                    all_contacts[particles_set] = canonical_contact

    return all_contacts
