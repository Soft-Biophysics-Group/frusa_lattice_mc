// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#ifndef MODEL_HEADER_H
#define MODEL_HEADER_H

#include <iostream>
#include <fstream>

/*Select the model library*/
#if defined DEFAULT
#include "default_include.h"
#elif defined FIELDS
#include "fields_parameters.h"
#include "fields_state.h"
#include "fields_geometry.h"
#include "fields_interactions.h"
#include "fields_update.h"
using namespace fields_space;
#endif

namespace model_space{
  /*Definition of model class*/
  class model {
    private:
      /*
       * Private variables
       */

      // Parameters from the input file
      model_parameters_struct parameters;

      // Structure containing the information about the current state
      // of the system
      state_struct state;

      // Structure containing the information about the interactions in the 
      // current state of the system
      interactions_struct interactions;

      // Structure containing information about MC averages
      averages_struct averages;

    public:

      /*Class constructor*/
      model();
 
      /*
       * Required public routines of the class
       */

      // Print the information about the current state of the system
      void print_model_state();

      // Save the current state of the system to a file "state_output"
      void save_model_state(std::string state_output);

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
       
  };
}
#endif
