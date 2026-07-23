"""Connected-component decompositions of a LatticeState."""

from dataclasses import dataclass
from functools import cached_property
from typing import Callable, Container, Sequence, TypeAlias

import numpy as np
from numpy.typing import NDArray

from ..lattice_state import LatticeState
from typing import Self

SitePair: TypeAlias = frozenset[int]
Contact: TypeAlias = tuple[int, int]
PairsToContacts: TypeAlias = dict[frozenset[int], Contact]


@dataclass(eq=False)
class Cluster:
    """One connected component. Built by `Clusters`, not directly."""

    sites: frozenset[int]
    state: LatticeState
    contacts: PairsToContacts

    @property
    def size(self) -> int:
        return len(self.sites)

    def count_contacts(self, contact_to_count: Container[Contact]) -> int:
        return sum(1 for c in self.contacts.values() if c in contact_to_count)

    @cached_property
    def n_interfaces_with_exterior(self) -> int:
        """Bonds from this cluster to an empty site."""
        full = self.state.full_sites_set
        return sum(
            1
            for site in self.sites
            for nb in self.state.lattice.get_neighbour_sites(site)
            if nb not in full
        )

    @cached_property
    def outer_perimeter(self) -> int:
        """Bonds from the cluster to the outside, not counting bonds towards holes."""
        latt = self.state.lattice

        # Seed in a cluster-free column: a hole cannot enclose a column, so it is outside
        occupied_cols = {latt.lattice_site_to_lattice_coords(s)[0] for s in self.sites}
        free_col = next((x for x in range(latt.lx) if x not in occupied_cols), None)
        if free_col is None:
            raise ValueError("cluster spans every column: it has no outside")
        outside_seed = latt.lattice_coords_to_lattice_site(free_col, 0, 0)

        # Flood fill the exterior. Holes are never reached, so their bonds are not counted
        outside = {outside_seed}
        stack = [outside_seed]
        while stack:
            site = stack.pop()
            for nb in latt.get_neighbour_sites(site):
                if nb not in self.sites and nb not in outside:
                    outside.add(nb)
                    stack.append(nb)

        return sum(
            1
            for site in self.sites
            for nb in latt.get_neighbour_sites(site)
            if nb in outside
        )

    @property
    def _relative_coords_cartesian(self) -> NDArray[np.float64]:
        r = np.array(list(self.coords_relative_to_center.values()))
        return r @ self.state.lattice.lattice_vectors_in_cartesian.T

    @property
    def bounding_box_dims(self) -> NDArray[np.float64]:
        return np.ptp(self._relative_coords_cartesian[:, :2], axis=0)

    @property
    def aspect_ratio_bb(self) -> float:
        """Measure the aspect ratio of a two-dimensional aggregate, taken as the ratio of the
        bounding box dimensions"""
        # Obtain the bounding box dimensions of the cluster
        bb_dims = self.bounding_box_dims
        max_dim = np.max(bb_dims)
        min_dim = np.min(bb_dims)
        return max_dim / min_dim

    @property
    def aspect_ratio_gyr(self) -> float:
        """Measure the aspect ratio of a two-dimensional aggregate, taken as the ratio between
        the eigenvalues of the gyration tensor"""
        latt = self.state.lattice
        # percolating clusters have no well-defined shape — guard first
        span = np.ptp(np.array(list(self.particle_coords_no_pbc.values())), axis=0)
        if (span[:2] >= np.array([latt.lx, latt.ly])).any():
            return np.nan

        r = self._relative_coords_cartesian
        G = (r.T @ r) / len(r)
        lam = np.linalg.eigvalsh(G)
        lam = lam[lam > 1e-9]  # drop flat dimensions (z in a 2D system)
        if lam.size < 2:
            return 1.0  # single point / line
        return float(np.sqrt(lam.max() / lam.min()))

    @property
    def particle_coords_no_pbc(self) -> dict[int, NDArray[np.int_]]:
        """Absolute (unwrapped) lattice coordinates for every site in the cluster,
        keyed by site index."""
        latt = self.state.lattice
        first = next(iter(self.sites))
        bond_sequences = {first: []}
        to_be_added = [first]

        for site in to_be_added:
            for i_bond, neigh in enumerate(latt.get_neighbour_sites(site)):
                if neigh not in self.sites or neigh in bond_sequences:
                    continue
                bond_sequences[neigh] = bond_sequences[site] + [i_bond]
                to_be_added.append(neigh)

        origin = latt.lattice_site_to_lattice_coords(first)
        return {
            site: origin
            + (
                np.sum([latt.bonds[i] for i in seq], axis=0)
                if seq
                else np.zeros(3, dtype=int)
            )
            for site, seq in bond_sequences.items()
        }

    @property
    def center_no_pbc(self) -> NDArray[np.int_]:
        raw_coords = np.vstack(
            [coord for coord in self.particle_coords_no_pbc.values()]
        )
        return np.mean(raw_coords, axis=0)

    @property
    def coords_relative_to_center(self) -> dict[int, NDArray[np.int_]]:
        """Lattice coordinates of each constitutive site relative to the cluster center"""
        return {
            site: coords - self.center_no_pbc
            for site, coords in self.particle_coords_no_pbc.items()
        }

def find_connected_site_sets(
    state: LatticeState,
    same_component: Callable[[int, int], bool],
    within: frozenset[int] | None = None,
) -> list[frozenset[int]]:
    """Site sets of the connected components of state, grouped by same_component.

    same_component is lambda x, y: True for aggregates, or an equality test on
    orientations for crystalline domains. If within is given, only those sites are
    considered, so the components returned are subsets of it.
    """
    full = state.full_sites_set if within is None else within
    # The tolist() converts the elements of full_sites from numpy int64 to int
    seeds = state.full_sites.tolist() if within is None else sorted(within)
    visited: set[int] = set()
    components: list[frozenset[int]] = []

    for seed in seeds:
        if seed in visited:
            continue
        # Build a cluster from an unvisited site
        component: set[int] = set()
        to_visit = {seed}
        visited.add(seed)
        while to_visit:
            site = to_visit.pop()
            component.add(site)
            for neigh in state.lattice.get_neighbour_sites(site):
                # Site gets added to component if if verifies same_component
                if (
                    neigh in full
                    and neigh not in visited
                    and same_component(site, neigh)
                ):
                    visited.add(neigh)
                    to_visit.add(neigh)
        components.append(frozenset(component))

    return components


class Clusters:
    """A decomposition of one LatticeState into connected components."""

    # Hidden attribute: the specific type of cluster I use I am planning on specializing
    # Cluster to different types (e.g. vortex assembly, pattern bulk...) with specific
    # properties. Having this _cluster_class property allows me to integrate these
    # particularities in the aggregate Clusters object
    _cluster_class: type[Cluster] = Cluster

    def __init__(self, state: LatticeState, site_sets: Sequence[frozenset[int]]):
        self.state = state
        self._site_sets = list(site_sets)

    @cached_property
    def clusters(self) -> list[Cluster]:
        internal, _ = self._partition
        return [
            self._cluster_class(sites, self.state, contacts)
            for sites, contacts in zip(self._site_sets, internal)
        ]

    @classmethod
    def connected_components_with_condition(cls, state, same_component) -> Self:
        return cls(state, find_connected_site_sets(state, same_component))

    # Let's specialize the connected_components to the 2 most common cases:
    @classmethod
    def connected_components(cls, state: LatticeState):
        """Builds all the connected components of state, i.e. the aggregates formed by connected
        particles with any orientation.
        """
        return cls.connected_components_with_condition(state, lambda x, y: True)

    # Let's specialize the connected_components to the 2 most common cases:
    @classmethod
    def crystalline_domains(cls, state: LatticeState):
        """Builds all the crystalline domains of state, i.e. the aggregates or parts of
        aggregates which have the same orientations.
        """
        orientations = state.orientations
        return cls.connected_components_with_condition(
            state, lambda x, y: orientations[x] == orientations[y]
        )

    def to_sub_clusters(self, same_component: Callable[[int, int], bool]) -> "Clusters":
        """Refines each cluster with the stricter condition same_component. Typical use
        case: break a patterned bulk into crystalline domains."""
        return Clusters(
            self.state,
            [
                sub_set
                for site_set in self._site_sets
                for sub_set in find_connected_site_sets(
                    self.state, same_component, within=site_set
                )
            ],
        )

    @cached_property
    def _partition(self) -> tuple[list[PairsToContacts], PairsToContacts]:
        """Assign every contact to its cluster, or to the boundary.
        Returns a list of Contacts, whose ith entry corresponds to the contacts within cluster
        i, and a Contacts instance listing the contacts at the boundaries between different
        clusters"""
        owner_of_site = {site: i for i, s in enumerate(self._site_sets) for site in s}
        internal: list[PairsToContacts] = [{} for _ in self._site_sets]
        boundary: PairsToContacts = {}
        for pair, contact in self.state.face_face_contacts.items():
            owners = {owner_of_site.get(site) for site in pair}
            # Both particles in the pair belong to the same component
            if len(owners) == 1 and None not in owners:
                # ty complains about this line because pop() returns None when the list is empty
                # Due to the "if" above, it's guaranteed to never be
                internal[owners.pop()][pair] = contact
            # Particles in the pair belong to different components
            else:
                boundary[pair] = contact
        return internal, boundary

    @cached_property
    def boundary_contacts(self) -> PairsToContacts:
        """Contacts between different clusters. Empty for aggregates, grain
        boundaries for crystalline domains."""
        return self._partition[1]

    def __len__(self) -> int:
        return len(self._site_sets)

    def __iter__(self):
        return iter(self.clusters)

    def __getitem__(self, i: int) -> Cluster:
        return self.clusters[i]

    @property
    def sizes(self) -> list[int]:
        return [len(s) for s in self._site_sets]

    def count_contacts_by_type(
        self, contact_types: Sequence[Sequence[tuple[int, int]]]
    ) -> NDArray[np.int_]:
        """(n_clusters, n_types) counts. Canonicalises each type class once.

        Given a set of contact categories in contact_types, counts the number of contacts of
        each category in each cluster.
        """
        canonical_contacts = [
            {self.state.particle.get_canonical_contact(*c) for c in group}
            for group in contact_types
        ]
        counts = np.zeros((len(self), len(canonical_contacts)), dtype=int)
        for i, cluster in enumerate(self):
            for j, types in enumerate(canonical_contacts):
                counts[i, j] = cluster.count_contacts(types)
        return counts


def get_aggregates(state: LatticeState) -> Clusters:
    return Clusters.connected_components(state)


def get_crystalline_domains(state: LatticeState) -> Clusters:
    return Clusters.crystalline_domains(state)


