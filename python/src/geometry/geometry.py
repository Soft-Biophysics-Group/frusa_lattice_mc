import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import config as cfg


class ParticleGeometry:
    def __init__(
        self,
        parameters_dict,
        lattice_spacing=1.0,
    ):
        self.permutations = parameters_dict["permutations"]
        self.opposite_bonds = parameters_dict["opposite_bonds"]
        self.bonds = parameters_dict["bonds"]
        self.n_faces = parameters_dict["n_faces"]
        self.lattice_vectors_in_cartesian = parameters_dict[
            "lattice_vectors_in_cartesian"
        ]
        self.ndims = parameters_dict["ndims"]
        self.rotation_matrix = parameters_dict["rotation_matrix"]
        self.lattice_spacing = lattice_spacing

    # ----- GEOMETRY -----
    def lattice_to_cartesian(
        self, x_lattice: int, y_lattice: int, z_lattice: int = 1
    ) -> int:
        if self.ndims == 2:
            return (
                x_lattice * self.lattice_vectors_in_cartesian[:, 0]
                + y_lattice * self.lattice_vectors_in_cartesian[:, 1]
            )
        if self.ndims == 3:
            return (
                x_lattice * self.lattice_vectors_in_cartesian[:, 0]
                + y_lattice * self.lattice_vectors_in_cartesian[:, 1]
                + z_lattice * self.lattice_vectors_in_cartesian[:, 2]
            )

    def rotate_point(self, coords, n_steps=1):
        """
        Applies `n_steps` rotation steps defined by `self.rotation_matrix` to cartesian
        coordinates coords in the clockwise direction.
        """
        for i in range(n_steps):
            coords = np.matmul(self.rotation_matrix, coords)
        return coords

    def get_bond(self, bond_vector):
        """
        If the `bond_vector`, a vector in lattice coordinates, links two neighbouring sites,
        this returns an unique index indicating the direction of this "bond"
        """
        matching_coeffs = self.bonds == bond_vector
        matching_vectors = np.prod(matching_coeffs, axis=1)
        if matching_vectors.sum() != 1:
            print("The vector does not correspond to a bond!")
            return -1
        else:
            return np.argmax(matching_vectors)

    def get_neighbour_sites_2d(self, site_index, lx, ly):
        """
        Return a 1D array of site indices corresponding to the neighbours of site_index.
        The jth element of this array is the neighbour of site_index following the jth bond.
        """
        x, y = lattice_site_to_lattice_coords_2d(site_index, lx)
        neighbours = []
        for bond in self.bonds:
            x_neighbour, y_neighbour = np.array([x, y]) + bond
            # Quick and dirty implementation of periodic boundary conditions
            if x_neighbour >= lx:
                x_neighbour -= lx
            if x_neighbour < 0:
                x_neighbour += lx
            if y_neighbour >= ly:
                y_neighbour -= ly
            if y_neighbour < 0:
                y_neighbour += ly

            neighbour_index = lattice_coords_to_lattice_site_2d(
                x_neighbour, y_neighbour, lx
            )
            neighbours.append(neighbour_index)
        return neighbours

    def get_neighbour_sites_3d(self, site_index, lx, ly, lz):
        x, y, z = lattice_site_to_lattice_coords_3d(site_index, lx, ly)
        neighbours = []
        for bond in self.bonds:
            x_neighbour, y_neighbour, z_neighbour = np.array([x, y, z]) + bond
            # Quick and dirty implementation of periodic boundary conditions
            if x_neighbour >= lx:
                x_neighbour -= lx
            if x_neighbour < 0:
                x_neighbour += lx
            if y_neighbour >= ly:
                y_neighbour -= ly
            if y_neighbour < 0:
                y_neighbour += ly
            if z_neighbour >= lz:
                z_neighbour -= lz
            if z_neighbour < 0:
                z_neighbour += lz

            neighbour_index = lattice_coords_to_lattice_site_2d(
                x_neighbour, y_neighbour, lx
            )
            neighbours.append(neighbour_index)
        return neighbours


# Helper functions in 2D
def lattice_site_to_lattice_coords_2d(site_index: int, lx: int):
    y_lattice = site_index // lx
    x_lattice = site_index - lx * y_lattice
    return np.array([x_lattice, y_lattice])


def lattice_coords_to_lattice_site_2d(x_lattice: int, y_lattice: int, lx: int):
    return x_lattice + y_lattice * lx

# Helper functions in 2D
def lattice_site_to_lattice_coords_3d(site_index: int, lx: int, ly:int):
    z_lattice = site_index // ly // lx
    y_lattice = ( site_index - lx * ly * z_lattice) // lx
    x_lattice = site_index - ly * lx * z_lattice - lx * y_lattice
    return np.array([x_lattice, y_lattice, z_lattice])


def lattice_coords_to_lattice_site_3d(
    x_lattice: int, y_lattice: int, z_lattice: int, lx: int, ly: int
):
    return x_lattice + y_lattice * lx + z_lattice * lx * ly


def get_full_sites_characteristics(results):
    """
    Taking in a results array which has the same format as the one returned by the c++ program,
    returns a 3D array with the occupied site as the first column, the particle type as the
    second, and the particle orientation as the third.
    """
    full_sites = []
    for site, (ptype, orientation) in enumerate(results.T):
        if orientation != -1:
            full_sites.append([site, ptype, orientation])
    return np.vstack(full_sites)

