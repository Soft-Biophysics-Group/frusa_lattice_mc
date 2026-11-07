"""Cluster.outer_perimeter on hand-built states with known answers: holes, pockets sealed
by two clusters, and sub-clusters touching other particles."""

import numpy as np
import pytest

from frusa_lattice_mc.analysis.clusters import Cluster
from frusa_lattice_mc.geometry.lattice_geometry import SquareLattice, TriangularLattice
from frusa_lattice_mc.geometry.particle_geometry import SquareParticle, TriangularParticle
from frusa_lattice_mc.lattice_state import LatticeState

L = 8


def make_state(lattice_cls, particle_cls, occupied_xy):
    lattice = lattice_cls(L, L)
    config = np.full((2, lattice.n_sites), -1, dtype=int)
    config[0, :] = 0
    for x, y in occupied_xy:
        config[1, lattice.lattice_coords_to_lattice_site(x, y)] = 0
    return LatticeState._from_lattice_config(config, lattice, particle_cls())


def make_cluster(state, cluster_xy):
    sites = frozenset(
        state.lattice.lattice_coords_to_lattice_site(x, y) for x, y in cluster_xy
    )
    return Cluster(sites, state, contacts={})


def square_cluster(cluster_xy, others_xy=()):
    state = make_state(SquareLattice, SquareParticle, [*cluster_xy, *others_xy])
    return make_cluster(state, cluster_xy)


RING_3X3 = [(x, y) for x in range(1, 4) for y in range(1, 4) if (x, y) != (2, 2)]


@pytest.mark.parametrize(
    "cluster_xy, expected",
    [
        ([(3, 3)], 4),
        ([(3, 3), (4, 3), (3, 4), (4, 4)], 8),
        # Straddles both periodic boundaries
        ([(7, 7), (0, 7), (7, 0), (0, 0)], 8),
    ],
)
def test_hole_free_matches_all_empty_interfaces(cluster_xy, expected):
    cluster = square_cluster(cluster_xy)
    assert cluster.outer_perimeter == expected
    assert cluster.n_interfaces_with_exterior == expected


def test_hole_excluded():
    cluster = square_cluster(RING_3X3)
    assert cluster.outer_perimeter == 12
    assert cluster.n_interfaces_with_exterior == 16


@pytest.mark.parametrize("include_other_clusters", [True, False])
def test_particle_inside_hole_excluded(include_other_clusters):
    cluster = square_cluster(RING_3X3, others_xy=[(2, 2)])
    assert cluster.count_outer_perimeter(include_other_clusters) == 12


@pytest.mark.parametrize("include_other_clusters", [True, False])
def test_pocket_sealed_by_two_clusters_is_outside(include_other_clusters):
    # Two L shapes touching only at corners enclose (2, 2) from the 4-neighbour fill
    cluster = square_cluster(
        [(1, 1), (2, 1), (1, 2)], others_xy=[(3, 3), (3, 2), (2, 3)]
    )
    assert cluster.count_outer_perimeter(include_other_clusters) == 8


def test_domain_buried_in_aggregate():
    others = [(x, y) for x in range(1, 4) for y in range(1, 4) if (x, y) != (2, 2)]
    cluster = square_cluster([(2, 2)], others_xy=others)
    assert cluster.outer_inner_perimeter == 4
    assert cluster.count_outer_perimeter(include_other_clusters=False) == 0


def test_domain_on_aggregate_edge():
    cluster = square_cluster([(1, 1), (1, 2)], others_xy=[(2, 1), (2, 2)])
    assert cluster.outer_inner_perimeter == 6
    assert cluster.count_outer_perimeter(include_other_clusters=False) == 4


def test_triangular_hexagonal_ring():
    state = make_state(TriangularLattice, TriangularParticle, [])
    latt = state.lattice
    center = latt.lattice_coords_to_lattice_site(3, 3)
    ring = latt.get_neighbour_sites(center)
    ring_xy = [tuple(latt.lattice_site_to_lattice_coords(s)[:2]) for s in ring]
    cluster = make_cluster(make_state(TriangularLattice, TriangularParticle, ring_xy), ring_xy)
    # Each ring site: 2 ring neighbours, the hole, 3 outside
    assert cluster.outer_perimeter == 18
    assert cluster.n_interfaces_with_exterior == 24


def test_percolating_cluster_raises():
    stripe = [(x, 3) for x in range(L)]
    with pytest.raises(ValueError, match="percolating"):
        square_cluster(stripe).outer_perimeter
