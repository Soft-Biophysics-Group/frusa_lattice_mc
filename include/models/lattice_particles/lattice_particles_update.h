#ifndef LATTICEPARTICLES_UPDATE_HEADER_H
#define LATTICEPARTICLES_UPDATE_HEADER_H

#include "json.hpp"
#include "lattice_particles_geometry.h"
#include "lattice_particles_interactions.h"
#include "lattice_particles_parameters.h"
#include "lattice_particles_state.h"

namespace lattice_particles_space {
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
void update_system(state_struct &state, interactions_struct &interactions,
                   model_parameters_struct &parameters, double T);


/*
 * End of the required definitions for the model class
 */

/*
 * Library-specific definitions
 */

// Pick one move at random
mc_moves pick_random_move(model_parameters_struct &parameters);

// Returns the index of a random lattice site containing a particle
std::size_t select_random_full_index(state_struct &state,
                                     model_parameters_struct &parameters);

// Returns the index of a random lattice site containing no particle
std::size_t select_random_empty_index(state_struct &state,
                                      model_parameters_struct &parameters);

/*
 * The following functions perform a MC move and return the corresponding
 * The following functions perform a MC move and return the correspondingenergy
 * difference.
 */
// TODO Finish writing this once I'm done with minor redesigns
double attempt_swap_sites(int index1, int index2, state_struct &state,
                          model_parameters_struct &parameters,
                          interactions_struct &interactions, double T);
double attempt_swap_empty_full(state_struct &state, model_parameters_struct &parameters,
                       interactions_struct &interactions, double T);
double attempt_swap_full_full(state_struct &state, model_parameters_struct &parameters,
                      interactions_struct &interactions, double T);
double attempt_rotate(state_struct &state, model_parameters_struct &parameters,
                      interactions_struct &interactions, double T);
double attempt_mutate(state_struct &state, model_parameters_struct &parameters,
                      interactions_struct &interactions, double T);

// Helper function to determine if sites are neighbours
// bool are_neighbours(site_state& site_1, site_state& site_2);
bool are_neighbours(int site_1_index, int site_2_index,
                    interactions_struct &interactions, state_struct &state);

double neighbour_correction(int site_1_index, int site_2_index,
                            interactions_struct &interactions,
                            state_struct &state,
                            model_parameters_struct &parameters);

// Accept or reject a move associated with energy delta_e at temperature T
// according to the Metropolis-Hastings rule
bool is_move_accepted(double delta_e, double T,
                      model_parameters_struct &parameters);
} // namespace lattice_particles_space

#endif
