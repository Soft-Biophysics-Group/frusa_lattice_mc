#include "triangular.h"
#include "geometry.h"

namespace geometry_space {
namespace triangular_space {
int get_neighbour_triangular(const int site_ind, const int bond_ind,
                             const int lx, const int ly) {

  int i{};
  int j{};
  int neigh_ind{0};

  array_space::r_to_ij(site_ind, i, j, lx, ly);

  int im = array_space::mod(i - 1, lx);
  int ip = array_space::mod(i + 1, lx);
  int jm = array_space::mod(j - 1, ly);
  int jp = array_space::mod(j + 1, ly);

  switch (bond_ind) {
  case 0:
    array_space::ij_to_r(neigh_ind, ip, j, lx, ly);
    break;
  case 1:
    array_space::ij_to_r(neigh_ind, ip, jp, lx, ly);
    break;
  case 2:
    array_space::ij_to_r(neigh_ind, i, jp, lx, ly);
    break;
  case 3:
    array_space::ij_to_r(neigh_ind, im, j, lx, ly);
    break;
  case 4:
    array_space::ij_to_r(neigh_ind, im, jm, lx, ly);
    break;
  case 5:
    array_space::ij_to_r(neigh_ind, i, jm, lx, ly);
    break;
  }

  return neigh_ind;
}
} // namespace triangular_space
} // namespace geometry_space
