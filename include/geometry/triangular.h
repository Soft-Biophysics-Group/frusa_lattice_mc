#ifndef GEOMETRY_TRIANGULAR_H
#define GEOMETRY_TRIANGULAR_H

namespace geometry_space {

namespace triangular_space {

class triangular_lattice {
public:
  triangular_lattice(int lx, int ly, int n_orientations);

private:
  int n_orientations_m{};
  int n_edges_m{6};
  int l_x_m{};
  int l_y_m{};
};

} // namespace triangular_space

} // namespace geometry_space

#endif
