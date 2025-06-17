"""Generic class for describing lattice geometry.

Typically not used on its own, but rather through its lattice-specific children classes.
"""
import numpy as np
from numpy.typing import ArrayLike, NDArray
import config as cfg
from collections.abc import Mapping, MutableMapping
from scipy.spatial.transform import Rotation as R
from pathlib import Path
from . import cubic, triangular, fcc
from .particle_geometry import ParticleGeometry

from typing import TypeAlias
Bond: TypeAlias = tuple[int, int, int]

class LatticeGeometry:
    """Generic class for lattice geometry representation.

    Lattices are made of sites with unique indices, which are the hashing of the site location
    in lattice coordinates.
    Each site is linked to its nearest neighbours with bond vectors, which are given in an
    arbitrary, but consistent, order.

    Attributes:
        lattice_vectors_in_cartesian: 1D array of floats. Lattice basis vectors in cartesian
            coordinates
        lattice_spacing: float. Distance between 2 nearest neighbours in Cartesian coordinates
        bonds: list of tuples of 3 ints. Element i is the i-th vector linking one lattice site
            with one of its nenearest neighbours, in lattice coordinates.
        bond_to_bond_index: dict mapping a bond (tuple of 3 ints) to its arbitrary index, as
            defined in bonds.
        lx: int. Number of lattice sites in the x direction in lattice coordinates.
        ly: int. Number of lattice sites in the y direction in lattice coordinates.
        lz: int. Number of lattice sites in the z direction in lattice coordinates.
    """

    # New way I'm trying out of specializing the class.
    # See: https://peps.python.org/pep-0487/#subclass-registration
    _lattices: MutableMapping[str, type["LatticeGeometry"]] = {}

    def __init__(
        self,
        basis_vectors: ArrayLike | None = None,
        bonds: list[Bond] | None = None,
        lx: int = 1,
        ly: int = 1,
        lz: int = 1,
        lattice_spacing: float = 1.0,
    ):
        """Generic constructor for LatticeGeometry class.

        Args:
            basis_vectors: array of floats. Basis vectors of the lattice.
            bonds: list of bonds (tuple of 3 ints). List of vectors linking one site to its
                nearest neighbours, in lattice coordinates.
            lx: dimensions of the lattice in the x direction, in lattice coordinates.
            ly: dimensions of the lattice in the y direction, in lattice coordinates.
            lz: dimensions of the lattice in the z direction, in lattice coordinates.
                Optional, defaults to 1 (for 2D lattices)
            lattice_spacing: optional float. Cartesian distance between 2 lattice nearest
                neighbours. Defaults to 1.0
        """
        self.lattice_vectors_in_cartesian: NDArray[np.float64] = np.array(
            basis_vectors, dtype=np.float64
        )
        self.lattice_spacing: float = lattice_spacing
        self.bonds: list[Bond]
        if bonds is not None:
            self.bonds = bonds
        else:
            self.bonds = [(1, 0, 0)]
        self.bond_to_bond_index: Mapping[Bond, int] = {
            bond: i for i, bond in enumerate(self.bonds)
        }
        self.lx: int = lx
        self.ly: int = ly
        self.lz: int = lz

        self.particle_geometry: ParticleGeometry
    # ----- SUBCLASSES -----
    # Magic code to register lattices upon creation:
    # inspired by https://stackoverflow.com/a/52433482
    @classmethod
    def get_implemented_lattices(cls):
        return list(cls._lattices)

    def __init_subclass__(cls, lattice: str) -> None:
        cls.lattice: str = lattice
        cls._lattices[lattice] = cls

    @classmethod
    def from_lattice_name(
        cls,
        lattice_name: str,
        lx: int,
        ly: int,
        lz: int = 1,
        lattice_spacing: float = 1.0,
    ):
        if lattice_name in cls._lattices.keys():
            inst =  cls._lattices[lattice_name](
                lx=lx, ly=ly, lz=lz, lattice_spacing=lattice_spacing
            )
            inst.particle_geometry = ParticleGeometry.from_lattice_name(lattice_name)
            return inst
        else:
            raise RuntimeError("Invalid lattice type selected")

    @classmethod
    def from_model_file(
        cls,
        model_file: str | Path = cfg.default_model_params_file,
        lattice_spacing: float = 1.0,
    ):
        model_params = cfg.load_model_file(model_file)
        lx: int = model_params["lx"]
        ly: int = model_params["ly"]
        lz: int = model_params["lz"]
        lattice_name = model_params["lattice_name"]

        return cls.from_lattice_name(lattice_name, lx, ly, lz, lattice_spacing)

    # ----- GEOMETRY -----
    def lattice_to_cartesian(
        self, x_lattice: int, y_lattice: int, z_lattice: int = 1
    ) -> NDArray[np.float_]:
        return (
            x_lattice * self.lattice_vectors_in_cartesian[:, 0]
            + y_lattice * self.lattice_vectors_in_cartesian[:, 1]
            + z_lattice * self.lattice_vectors_in_cartesian[:, 2]
        )

    def get_bond(self, bond: ArrayLike) -> int:
        """If the `bond_vector`, a vector in lattice coordinates, links two neighbouring sites,
        this returns an unique index indicating the direction of this "bond"
        """
        # Turn contents of bond into ints
        bond_arr = np.array(bond, dtype=int)
        # print(f"before pbc {bond_arr}")
        bond_arr = self.apply_pbc_vector(bond_arr)
        # print(f"after pbc {bond_arr}")
        # And then into tuple for hashing
        bond_tup = tuple(bond_arr)
        if bond_tup not in self.bond_to_bond_index.keys():
            return -1
        else:
            return self.bond_to_bond_index[bond_tup]

    def get_bond_from_sites(self, site_1: int, site_2: int) -> int:
        """Attempts to find the bond index between site_1 and site_2, taking site_1 as the
        center of frame of reference. Returns -1 if the two sites are not neighbours.
        """
        site_1_lattice_coords = self.lattice_site_to_lattice_coords(site_1)
        site_2_lattice_coords = self.lattice_site_to_lattice_coords(site_2)
        bond_coords = site_2_lattice_coords - site_1_lattice_coords
        print(bond_coords)

        # bond_coords = np.array(self.apply_pbc(*bond_coords))
        print(bond_coords)
        bond = self.get_bond(bond_coords)

        return bond

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
        """Return as 1D array of site indices corresponding to the neighbours of site_index.
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
            x_neighbour, y_neighbour, z_neighbour = self.apply_pbc(
                x_neighbour, y_neighbour, z_neighbour
            )
            neighbour_index = self.lattice_coords_to_lattice_site(
                x_neighbour, y_neighbour, z_neighbour
            )
            neighbours.append(neighbour_index)
        return neighbours

    def apply_pbc(self, x: int, y: int, z: int) -> tuple[int, int, int]:
        x = np.mod(x, self.lx)
        y = np.mod(y, self.ly)
        z = np.mod(z, self.lz)

        return x, y, z

    def apply_pbc_vector(self, vec: NDArray[np.float_]):
        """Applies periodic boundary conditions to a vector linking 2 sites.

        If the vector has components with magnitude larger than half of the lattice, we wrap
        these components back
        """

        vec_pbc = np.copy(vec)
        x, y, z = vec_pbc
        if np.abs(x) >= self.lx // 2:
            vec_pbc[0] = -np.sign(x) * (self.lx - np.abs(x))
        if np.abs(y) >= self.ly // 2:
            vec_pbc[1] = -np.sign(y) * (self.ly - np.abs(y))
        if np.abs(z) >= self.lz // 2:
            vec_pbc[2] = -np.sign(z) * (self.lz - np.abs(z))

        return vec_pbc

    def get_faces_in_contact_and_bond(
        self,
        site_1_index: int,
        orientation_1: int,
        site_2_index: int,
        orientation_2: int,
    ):
        site_1_lattice_coords = np.array(
            self.lattice_site_to_lattice_coords(site_1_index), dtype=int
        )
        site_2_lattice_coords = np.array(
            self.lattice_site_to_lattice_coords(site_2_index), dtype=int
        )

        bond_lattice_coords = self.apply_pbc_vector(
            site_2_lattice_coords - site_1_lattice_coords
        )
        bond_index = self.get_bond(bond_lattice_coords)

        if bond_index == -1:
            return (-1, -1, -1)

        face_1, face_2 = self.particle_geometry.get_faces_in_contact(
            orientation_1, orientation_2, bond_index
        )

        return face_1, face_2, bond_index


# ----- LATTICE SPECIALIZATIONS -----


# The lz = 1.0 is here for conssistency. I put it last so that it doesn't get in the way.
class TriangularLattice(LatticeGeometry, lattice="triangular"):
    def __init__(
        self, lx: int = 1, ly: int = 1, lattice_spacing: float = 1.0, lz: int = 1
    ):
        super().__init__(
            triangular.BASIS_VECTORS, triangular.BONDS, lx, ly, 1, lattice_spacing
        )


class CubicLattice(LatticeGeometry, lattice="cubic"):
    def __init__(
        self, lx: int = 1, ly: int = 1, lz: int = 1, lattice_spacing: float = 1.0
    ):
        super().__init__(cubic.BASIS_VECTORS, cubic.BONDS, lx, ly, lz, lattice_spacing)


class FccLattice(LatticeGeometry, lattice="fcc"):
    def __init__(
        self, lx: int = 1, ly: int = 1, lz: int = 1, lattice_spacing: float = 1.0
    ):
        super().__init__(fcc.BASIS_VECTORS, fcc.BONDS, lx, ly, lz, lattice_spacing)
