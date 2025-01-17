#ifndef PARTICLES_UPDATE_HEADER_H
#define PARTICLES_UPDATE_HEADER_H

/**
 * Functions to update the state of particles occupying the sites of a Bravais
 * lattice according to the Metropolis-Hastings algorithm, and update the
 * energy of the system accordingly.
 */

#include "json.hpp"
#include "geometry.h"
#include "particles_interactions.h"
#include "particles_parameters.h"
#include "particles_state.h"

namespace particles_space {
using json = nlohmann::json;
/*
 * Definitions required for the public routines of the model class
 */

// Function used to perform state.N updates of the system (considered a
// single MC step)
// state        - configuration of the system before the update
// interactions - energetics of the system before the update
// parameters   - used to access random number generator (parameters.rng)
// T            - annealing temperature (not the same as T_model!)
void update_system(state_struct& state,
                   interactions_struct& interactions,
                   model_parameters_struct& parameters,
                   geometry_space::Geometry geometry,
                   double T);

/*
 * End of the required definitions for the model class
 */

/*
 * Library-specific definitions
 */

mc_moves pick_random_move(model_parameters_struct &parameters);

std::size_t select_random_full_index(state_struct &state,
                                     model_parameters_struct &parameters);
std::size_t select_random_empty_index(state_struct &state,
                                      model_parameters_struct &parameters);

// Function specifically to perform random rotations to avoid code reuse
// Returns the orientation of the particle before rotation
int perform_random_rotation(state_struct &state,
                            model_parameters_struct &parameters,
                            int site_index);

// Measure the total energy of a pair of particles, avoiding double counting if
// they occupy neighbouring sites
double measure_pair_energy(int index1,
                           int index2,
                           int bond,
                           state_struct& state,
                           interactions_struct& interactions,
                           geometry_space::Geometry& geometry);

/*
 * The following functions perform a MC move and return the corresponding energy
 * difference.
 */
double attempt_swap_sites(int index1,
                          int index2,
                          state_struct& state,
                          model_parameters_struct& parameters,
                          interactions_struct& interactions,
                          geometry_space::Geometry geometry,
                          double T);
double attempt_swap_empty_full(state_struct& state,
                               model_parameters_struct& parameters,
                               interactions_struct& interactions,
                               geometry_space::Geometry geometry,
                               double T);
double attempt_swap_full_full(state_struct& state,
                              model_parameters_struct& parameters,
                              interactions_struct& interactions,
                              geometry_space::Geometry geometry,
                              double T);
double attempt_rotate(state_struct& state,
                      model_parameters_struct& parameters,
                      interactions_struct& interactions,
                      geometry_space::Geometry geometry,
                      double T);
double attempt_mutate(state_struct& state,
                      model_parameters_struct& parameters,
                      interactions_struct& interactions,
                      geometry_space::Geometry geometry,
                      double T);
double attempt_rotate_and_swap_w_empty(state_struct& state,
                                       model_parameters_struct& parameters,
                                       interactions_struct& interactions,
                                       geometry_space::Geometry geometry,
                                       double T);

// Accept or reject a move associated with energy delta_e at temperature T
// according to the Metropolis-Hastings rule
bool is_move_accepted(double delta_e, double T,
                      model_parameters_struct &parameters);
} // namespace lattice_particles_space

#endif
