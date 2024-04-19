#include "geometry.h"
#include "triangular.h"
#include <stdexcept>

namespace geometry_space {
Geometry::Geometry(lattice_options lattice, int lx, int ly, int lz)
    : lattice_m{lattice}, lx_m{lx}, ly_m{ly}, lz_m{lz} {
  switch (lattice_m) {
  case lattice_options::triangular:
    n_neighbours = 6;
    break;
  default:
    throw(std::runtime_error("Invalid lattice option"));
  }
}

int Geometry::get_neighbour(const int site_ind, const int bond_ind) const {
  switch (lattice_m) {
  case lattice_options::triangular:
    return triangular_space::get_neighbour_triangular(site_ind, bond_ind, lx_m,
                                                      ly_m);
  default:
    throw(std::runtime_error("Invalid lattice option"));
  }
}
} // namespace geometry_space
