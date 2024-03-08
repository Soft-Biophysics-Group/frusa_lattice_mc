#ifndef FIELDS_CHAIN_HEADER_H
#define FIELDS_CHAIN_HEADER_H

#include "fields_state.h"
#include "vector_utils.h"

namespace fields_space{
  /*
   * Routines used to initialize the couplings for a 1d system of particles
   * with 2 orientations
   */
  vec3d get_coupling_matrix(vec1d couplings);

  vec1i get_neighbours(int r, state_struct &state);
}

#endif
