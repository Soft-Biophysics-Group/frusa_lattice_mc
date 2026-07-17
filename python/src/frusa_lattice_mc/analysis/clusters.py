"""Connected-component decompositions of a LatticeState."""

from dataclasses import dataclass
from functools import cached_property
from typing import Callable, Container, Sequence, TypeAlias

import numpy as np
from numpy.typing import NDArray

from ..lattice_state import LatticeState

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


class Clusters:
    """A decomposition of one LatticeState into connected components."""

    def __init__(self, state: LatticeState, site_sets: Sequence[frozenset[int]]):
        self.state = state
        self._site_sets = list(site_sets)

    @classmethod
    def connected_components(
        cls,
        state: LatticeState,
        same_component: Callable[[int, int], bool] = lambda x, y: True,
    ) -> "Clusters":
        """Finds the connected components of a LatticeState sharing some conditions, same_component, on their orientation.
        same_condition is typically going to be lambda x,y : True for connected components,
        or lambda x, y: x==y for crystalline domains"""
        full = state.full_sites_set
        visited: set[int] = set()
        components: list[frozenset[int]] = []

        # The tolist() converts the elements of full_sites from numpy int64 to int
        for seed in state.full_sites.tolist():
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

        return cls(state, components)

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
    def clusters(self) -> list[Cluster]:
        internal, _ = self._partition
        return [
            Cluster(sites, self.state, contacts)
            for sites, contacts in zip(self._site_sets, internal)
        ]

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
    return Clusters.connected_components(state, lambda a, b: True)


def get_crystalline_domains(state: LatticeState) -> Clusters:
    orientations = state.orientations
    return Clusters.connected_components(
        state, lambda a, b: orientations[a] == orientations[b]
    )
