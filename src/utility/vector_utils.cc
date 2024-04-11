#include "vector_utils.h"

namespace array_space{

  int hash_into_state(int type, int orientation, int n_orientations) {
    if (orientation == 0)
      return 0;
    else
      return type * n_orientations + orientation;
  }
  int hash_3integers(int x1, int x2, int y, int x_range) {
      return x1 * x_range + x2 * x_range^2 + y;
  }
  int hash_into_contact(int state1, int state2, int edge1, int edge2,
                        int n_states) {
    // Here we need to determine the ordering of the 2 particles making the
    // pair for hashing
    // If one of the sites is empty, i.e. state 0, then it has to be the second
    // particle
    if (state1 == 0)
      return hash_3integers(state2, state1, edge2, n_states);
    else if (state2 == 0)
      return hash_3integers(state1, state2, edge1, n_states);
    // If both sides full, then the one with the lowest edge number will be
    // particle 1, and its edge is taken as reference
    else if (edge1 >= edge2)
      return hash_3integers(state1, state2, edge1, n_states);
    else
      return hash_3integers(state2, state1, edge2, n_states);
  }

  void ij_to_r(int &r, int i, int j, int Lx){
    r = i + Lx*j;
  }

  void r_to_ij(int r, int &i, int &j, int Lx){
    j = r/Lx;
    i = r-Lx*j;
  }

  void ijk_to_r(int &r, int i, int j, int k, int Lx, int Ly){
    r = i + Lx*j + Lx*Ly*k;
  }

  void r_to_ijk(int r, int &i, int &j, int &k, int Lx, int Ly){
    k = r/(Lx*Ly);
    j = (r-Lx*Ly*k)/Lx;
    i = r-Lx*j-Lx*Ly*k;
  }

  int mod(int a, int b){
    return (a % b + b) % b;
  }
}
