// Copyright (c) 2024 Andrey Zelenskiy
// Part of frusa_mc, released under BSD 3-Clause License.

#include "model.h"

namespace model_space{
  model::model(){
    /*
     * Initialize the system of particles and calculate the initial energy
     */

    initialize_state(state,parameters);
    
    initialize_interactions(state,interactions,parameters);
    
  }

  void model::print_model_state(){
    print_state(state);
  }

  void model::save_model_state(std::string state_output){
    save_state(state, state_output);
  }
  
  void model::print_model_interactions(){
    print_interactions(state,interactions);
  }

  void model::print_model_energy(){
    print_energy(state,interactions);
  }

  void model::update_model_system(double T){
    update_system(state,interactions,parameters,T);
  }

  void model::initialize_model_averages(){
    initialize_averages(averages,parameters);
  }

  void model::update_model_averages(double T){
    update_averages(averages,state,interactions,parameters,T);
  }

  void model::save_model_averages(double T, int mcs_av){
    save_averages(averages,state,parameters,T,mcs_av);
  }
}
