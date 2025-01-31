// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef MODEL_HEADER_H
#define MODEL_HEADER_H

#include <fstream>
#include <iostream>
#include <string>

/*Select the model library*/
// #include "default_include.h"

// #include "fields_parameters.h"
// #include "fields_state.h"
// #include "fields_geometry.h"
// #include "fields_interactions.h"
// #include "fields_update.h"

#include "geometry.h"
#include "particles_averages.h"
#include "particles_interactions.h"
#include "particles_parameters.h"
#include "particles_state.h"
#include "particles_update.h"
#include "particles_records.h"

namespace model_space {

enum model_options {
  fields,
  particles,
  n_models
};

static const inline std::array<std::string, model_options::n_models>
    lattice_str_arr{"fields", "particles"};

struct model_parameters_struct{
  model_parameters_struct(model_options model);
};

/*Definition of model class*/
class model {
private:
  /*
   * Private variables
   */

  // Geometry of the lattice
  geometry_space::Geometry geometry;

  // Parameters from the input file
  particles_space::model_parameters_struct parameters;

  // Structure containing the information about the current state
  // of the system
  particles_space::state_struct state;

  // Structure containing the information about the interactions in the
  // current state of the system
  particles_space::interactions_struct interactions;

  // Structure containing information about MC averages
  particles_space::averages_struct averages;

  // Structure containing records of energy after each lattice update
  particles_space::records_struct records;

public:
  /*Class constructor*/
  model(std::string& model_params_file);

  /*
   * Required public routines of the class
   */

  // Print the information about the current state of the system
  void print_model_state();

  // Save the current state of the system to a file "state_output"
  void save_model_state(std::string& state_output);

  // Print the information about interactions in the current state of the
  // system
  void print_model_interactions();

  void print_model_energy();

  // Update the state of the system at annealing temperature T
  void update_model_system(double T);

  // Initialize the containers to store the selected averages
  void initialize_model_averages();

  // Update the selected simulation averages
  void update_model_averages(double T);

  // Save the selected simulation averages to the corresponding files
  void save_model_averages(double T, int mcs_av);

  // Save a set of recorded energies after a lattice update
  void update_model_records();

  // Save a set of recorded energies after a lattice update
  void save_model_records(double T);
};
} // namespace model_space
#endif
