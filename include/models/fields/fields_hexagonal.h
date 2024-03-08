#ifndef FIELDS_HEXAGONAL_HEADER_H
#define FIELDS_HEXAGONAL_HEADER_H

#include "fields_state.h"
#include "fields_interactions.h"
#include "fields_update.h"

#include "vector_utils.h"

namespace fields_space{
  /*
   * Routines used to initialize the couplings for a 2d system of particles
   * with 6 orientations on a hexagonal lattice
   */
  vec3d get_coupling_matrix(vec1d couplings);

  vec1i get_neighbours(int r, state_struct &state);
}

#endif
