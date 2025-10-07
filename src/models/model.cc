// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#include "model.h"

namespace model_space
{
model_parameters_struct::model_parameters_struct(model_options model)
{
  switch (model) {
    case fields:;
    default:;
  }
}
model::model(std::string& model_params_file)
    : geometry {geometry_space::Geometry(model_params_file)}
    , parameters {particles_space::model_parameters_struct(model_params_file)}
    , state {particles_space::state_struct {}}
    , interactions {particles_space::interactions_struct {}}
    , averages {particles_space::averages_struct {}}
    , records {particles_space::records_struct {}}
{
  /*
   * Initialize the system of particles and calculate the initial energy
   */

  std::cout << "Got here" ;
  particles_space::initialize_state(state, parameters, geometry);

  particles_space::initialize_interactions(
      state, interactions, parameters, geometry);
}

void model::print_model_state()
{
  particles_space::print_state(state);
}

void model::save_model_state(std::string& state_output)
{
  particles_space::save_state(state, state_output);
}

void model::print_model_interactions()
{
  particles_space::print_interactions(interactions);
}

void model::print_model_energy()
{
  particles_space::print_energy(interactions);
}

void model::update_model_system(double T)
{
  particles_space::update_system(state, interactions, parameters, geometry, T);
}

void model::initialize_model_averages()
{
  initialize_averages(averages, parameters);
}

void model::update_model_averages(double T)
{
  update_averages(averages, state, interactions, parameters, T);
}

void model::save_model_averages(double T, int mcs_av)
{
  save_averages(averages, state, parameters, T, mcs_av);
}

void model::update_model_records()
{
  update_records(parameters, interactions, records);
}

void model::save_model_records(double T)
{
  save_records(parameters, T, records);
}

}  // namespace model_space
