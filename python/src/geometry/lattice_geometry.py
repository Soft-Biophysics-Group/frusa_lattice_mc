import numpy as np
from numpy.typing import ArrayLike, NDArray
import config as cfg

from typing import TypeAlias
Bond: TypeAlias = NDArray[np.int_]


class LatticeGeometry:
    def __init__(
        self,
        basis_vectors: ArrayLike,
        bonds: ArrayLike,
        lx: int,
        ly: int,
        lz: int = 1,
        lattice_spacing: float = 1.0,
    ):
        self.lattice_vectors_in_cartesian: NDArray[np.float64] = np.array(
            basis_vectors, dtype=np.float64
        )
        self.lattice_spacing: float = lattice_spacing
        self.bonds: NDArray[np.int_] = np.array(bonds, dtype = np.int_)
        self.lx: int = lx
        self.ly: int = ly
        self.lz: int = lz

    # ----- GEOMETRY -----
    def lattice_to_cartesian(
        self, x_lattice: int, y_lattice: int, z_lattice: int = 1
    ) -> ArrayLike:
        return (
            x_lattice * self.lattice_vectors_in_cartesian[:, 0]
            + y_lattice * self.lattice_vectors_in_cartesian[:, 1]
            + z_lattice * self.lattice_vectors_in_cartesian[:, 2]
        )

    def get_bond(self, bond: ArrayLike) -> int:
        """
        If the `bond_vector`, a vector in lattice coordinates, links two neighbouring sites,
        this returns an unique index indicating the direction of this "bond"
        """
        bond_vector: NDArray[np.int_] = np.array(bond, dtype=int)
        matching_coeffs: NDArray[np.bool_] = (self.bonds == bond_vector)
        matching_vectors: NDArray[np.bool_] = np.prod(matching_coeffs, axis=1)
        if matching_vectors.sum() != 1:
            print("The vector does not correspond to a bond!")
            return -1
        else:
            return int(np.argmax(matching_vectors))

    def lattice_site_to_lattice_coords(self, site_index: int) -> NDArray[np.int_]:
        z_lattice = site_index // self.ly // self.lx
        y_lattice = (site_index - self.lx * self.ly * z_lattice) // self.lx
        x_lattice = site_index - self.ly * self.lx * z_lattice - self.lx * y_lattice
        return np.array([x_lattice, y_lattice, z_lattice])

    def lattice_coords_to_lattice_site(
        self, x_lattice: int, y_lattice: int, z_lattice: int = 1
    ) -> int:
        return x_lattice + y_lattice * self.lx + z_lattice * self.lx * self.ly

    def get_neighbour_sites(self, site_index: int) -> list[int]:
        """
        Return a 1D array of site indices corresponding to the neighbours of site_index.
        The jth element of this array is the neighbour of site_index following the jth bond.
        """
        x: int
        y: int
        z: int
        x, y, z = self.lattice_site_to_lattice_coords(site_index)

        neighbours: list[int] = []
        bond: Bond
        for bond in self.bonds:
            x_neighbour: int
            y_neighbour: int
            z_neighbour: int
            x_neighbour, y_neighbour, z_neighbour = np.array([x, y, z]) + bond
            # Quick and dirty implementation of periodic boundary conditions
            if x_neighbour >= self.lx:
                x_neighbour -= self.lx
            if x_neighbour < 0:
                x_neighbour += self.lx
            if y_neighbour >= self.ly:
                y_neighbour -= self.ly
            if y_neighbour < 0:
                y_neighbour += self.ly
            if z_neighbour >= self.lz:
                z_neighbour -= self.lz
            if z_neighbour < 0:
                z_neighbour += self.lz

            neighbour_index = self.lattice_coords_to_lattice_site(
                x_neighbour, y_neighbour, z_neighbour
            )
            neighbours.append(neighbour_index)
        return neighbours


# def get_full_sites_characteristics(results):
#     """
#     Taking in a results array which has the same format as the one returned by the c++ program,
#     returns a 3D array with the occupied site as the first column, the particle type as the
#     second, and the particle orientation as the third.
#     """
#     full_sites = []
#     for site, (ptype, orientation) in enumerate(results.T):
#         if orientation != -1:
#             full_sites.append([site, ptype, orientation])
#     return np.vstack(full_sites)
