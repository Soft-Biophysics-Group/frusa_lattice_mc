#ifndef PARTICLES_RECORDS_H
#define PARTICLES_RECORDS_H

#include <fstream>
#include <iostream>
#include <string>

#include "particles_interactions.h"
#include "particles_parameters.h"
#include "particles_state.h"
#include "vector_utils.h"
#include "io_utils.h"

namespace particles_space
{
/*
 * Definitions required for the public routines of the model class
 */

// Structure containing the containers for the average quantities:
// e_record - energy after each lattice update
// e_record - temperature after each lattice update. Can be re-derived, but more
// convenient like this
struct records_struct
{
  vec1d e_records {};
};

void update_records(model_parameters_struct& parameters,
                    interactions_struct& interactions,
                    records_struct& records);

void save_records(model_parameters_struct& parameters,
                  double T,
                  records_struct& records);

/*
 * End of the required definitions for the model class
 */

}  // namespace particles_space

#endif
