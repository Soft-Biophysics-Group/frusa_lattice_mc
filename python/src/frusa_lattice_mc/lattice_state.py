"""Vincent Ouazan-Reboul, 2026-07-17
Unified utilities to load and manipulate lattice states. Useful for plotting, data analysis
"""

from pathlib import Path
import numpy as np
from numpy.typing import NDArray
from functools import cached_property
from .geometry import LatticeGeometry
from .geometry import ParticleGeometry


class LatticeState:
    """One simulation snapshot: the particle configuration plus the geometry it lives on.

    `lattice_config` is a (2, n_sites) array, row 0 being particle type and row 1 particle
    orientation. Orientation -1 marks an empty site, whatever the type.
    """

    lattice_config: NDArray[np.int_]
    lattice: LatticeGeometry
    particle: ParticleGeometry

    def __init__(self, path_to_config: str | Path, model_file: str | Path) -> None:
        self.lattice_config = np.loadtxt(path_to_config, dtype=int)
        self.lattice = LatticeGeometry.from_model_file(model_file)
        self.particle = ParticleGeometry.from_model_file(model_file)

    @classmethod
    def _from_lattice_config(
        cls,
        lattice_config: NDArray[np.int_],
        lattice: LatticeGeometry,
        particle: ParticleGeometry,
    ) -> "LatticeState":
        """Build a state from an in-memory configuration, reusing existing geometry.

        The geometry objects are immutable, so sharing them between states is safe and
        avoids reloading the model file.
        """
        new_state = cls.__new__(cls)
        new_state.lattice_config = lattice_config
        new_state.lattice = lattice
        new_state.particle = particle
        return new_state

    @classmethod
    def from_dir(
        cls,
        struct_folder: str | Path,
        model_file: str | Path,
        struct_index: int | None = None,
    ):
        """Pull a structure directly from a folder. If struct_index is not specified, pulls the
        final structure of the run"""
        struct_path = Path(struct_folder)
        if struct_index is not None:
            file_path = struct_path / f"structure_{struct_index}.dat"
        else:
            file_path = struct_path / "final_structure.dat"

        return cls(file_path, model_file)

    @cached_property
    def full_sites(self) -> NDArray[np.int64]:
        return np.where(self.lattice_config[1, :] != -1)[0]

    @cached_property
    def full_sites_set(self) -> set[int]:
        """Full sites as a set, for O(1) membership tests in the site loops."""
        return set(self.full_sites.tolist())

    @cached_property
    def orientations(self) -> NDArray[np.int_]:
        """Orientation of each site, indexed by site number. -1 marks an empty site."""
        return self.lattice_config[1, :]

    @cached_property
    def get_full_sites_characteristics(self):
        sites = self.full_sites
        types = self.lattice_config[0, sites]
        orientations = self.lattice_config[1, sites]
        return np.vstack([sites, types, orientations]).T

    @cached_property
    def face_face_contacts(self) -> dict[frozenset[int], tuple[int, int]]:
        """All the particle face-face contacts in this state, as a dict for data analysis.

        Returns:
        A dict mapping each pair of neighbouring full sites (as a frozenset of the two site
        indices) to the canonical form of the face-face contact between the two particles
        occupying them. Only pairs of full sites appear as keys.
        """

        all_contacts: dict[frozenset[int], tuple[int, int]] = {}

        for site_1 in self.full_sites:
            orientation_1 = self.orientations[site_1]
            neighbours_of_1 = self.lattice.get_neighbour_sites(site_1)
            for site_2 in neighbours_of_1:
                if site_2 in self.full_sites_set:
                    particles_set = frozenset((site_1, site_2))
                    if particles_set not in all_contacts:
                        orientation_2 = self.orientations[site_2]
                        face_1, face_2, _ = self.lattice.get_faces_in_contact_and_bond(
                            site_1, orientation_1, site_2, orientation_2
                        )
                        canonical_contact = self.particle.get_canonical_contact(
                            face_1, face_2
                        )
                        all_contacts[particles_set] = canonical_contact

        return all_contacts

    def _get_translated_site_index(
        self, site: int, translation_vec: list[int] | NDArray[np.int_]
    ) -> int:
        site_coords_lattice = self.lattice.lattice_site_to_lattice_coords(site)
        translated_coords = site_coords_lattice + translation_vec
        new_site_coords_lattice = self.lattice.apply_pbc(*translated_coords)
        new_site = self.lattice.lattice_coords_to_lattice_site(*new_site_coords_lattice)

        return new_site

    def translate(
        self,
        translation_vec: list[int] | NDArray[np.int_],
    ) -> "LatticeState":
        """Return a new state with every particle shifted by `translation_vec`.

        The current state is left untouched. The new state shares this one's geometry.
        """
        new_config = np.zeros_like(self.lattice_config) - 1

        for site in self.full_sites:
            new_site = self._get_translated_site_index(site, translation_vec)
            new_config[:, new_site] = self.lattice_config[:, site]

        return self._from_lattice_config(new_config, self.lattice, self.particle)
