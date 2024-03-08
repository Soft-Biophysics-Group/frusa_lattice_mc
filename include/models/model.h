#ifndef MODEL_HEADER_H
#define MODEL_HEADER_H

#include <iostream>
#include <fstream>

/*Select the model library*/
//#if defined FIELDS_1D_OPTION
#include "fields_chain.h"
using namespace fields_space;
//#else
//#include "default_model.h"
//#endif

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

      /*
       * Private routines of the class
       */

      // Function to initialize the system of anisotropic particles
      void initialize();
    
    public:

      /*Class constructor*/
      model();
 
      /*
       * Required public routines of the class
       */

      // Print the information about the current state of the system
      void print_model_state();

      // Save the current state of the system to a file
      void save_model_state();

      // Print the information about interactions in the current state of the 
      // system
      void print_model_interactions();

      // Update the state of the system
      void update_model_state(double);

      // Initialize the containers to store the selected averages
      void initialize_model_averages();

      // Update the selected simulation averages
      void update_model_averages(double);

      // Save the selected simulation averages to the corresponding files
      void save_model_averages();
       
  };
}
#endif
