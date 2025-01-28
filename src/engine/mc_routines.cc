// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#include "mc_routines.h"

namespace simulation_space{

  mc_parameters_struct::mc_parameters_struct(std::string_view mc_input){
    /*
     * Populate the struct using the input JSON file
     */

    std::ifstream mc_f;
    mc_f.open(mc_input);
    if(!mc_f){
      std::cerr << "Could not open JSON MC parameters file" << std::endl;
      exit(1);
    }

    json json_mc_params = json::parse(mc_f);

    mcs_eq            = json_mc_params["mcs_eq"].template get<int>();
    mcs_av            = json_mc_params["mcs_av"].template get<int>();
    cooling_schedule  =\
      json_mc_params["cooling_schedule"].template get<std::string>();
    Ti                = json_mc_params["Ti"].template get<double>();
    Tf                = json_mc_params["Tf"].template get<double>();
    Nt                = json_mc_params["Nt"].template get<int>();
    checkpoint_option =\
      json_mc_params["checkpoint_option"].template get<bool>();

    if(checkpoint_option){
      checkpoint_address =\
        json_mc_params["checkpoint_address"].template get<std::string>();
    }
    final_structure_address =\
      json_mc_params["final_structure_address"].template get<std::string>();
  }

  mc::mc(std::string_view mc_input): parameters {mc_parameters_struct(mc_input)}{

    // Map the cooling option on the integer variable
    try{
      if(parameters.cooling_schedule=="exponential"){
        cooling_option = 0;
      }
      else if(parameters.cooling_schedule=="linear"){
        cooling_option = 1;
      }
      else{
        throw parameters.cooling_schedule;
      }
    }

    catch(std::string cooling_schedule){
      std::cout << cooling_schedule << ": Incorrect cooling option!\n";
      exit(1);
    }

    // Define the array of temperatures for the annealing
    double dT = (parameters.Tf-parameters.Ti)/parameters.Nt;

    for(int i=0;i<parameters.Nt;i++){
      T_array.push_back(parameters.Ti+i*dT);
    }
  }

  void mc::print_mc_parameters(){
    std::cout << "\n---------------------------------------------\n";
    std::cout << "      Monte-Carlo simulation parameters\n";
    std::cout << "---------------------------------------------\n\n";

    std::cout << "Number of equilibration steps mcs_eq = ";
    std::cout << parameters.mcs_eq << "\n";
    std::cout << "Number of averaging steps mcs_av = ";
    std::cout << parameters.mcs_av << "\n\n";

    std::cout << "Initial temperature Ti = ";
    std::cout << parameters.Ti << "\n";
    std::cout << "Final temperature Tf = ";
    std::cout << parameters.Tf << "\n";
    std::cout << "Number of temperature steps Nt = ";
    std::cout << parameters.Nt << "\n\n";

    std::cout << "The selected cooling schedule is ";
    std::cout << parameters.cooling_schedule << "\n\n";

    if(parameters.checkpoint_option){
      std::cout << "Output location for checkpoints: ";
      std::cout << parameters.checkpoint_address << "\n\n";
    }

    std::cout << "Output location for the final structure: ";
    std::cout << parameters.final_structure_address << "\n\n";
  }

  void mc::t_scan(model_space::model &simulation_model){

    for (std::size_t i = 0; i < static_cast<std::size_t>(parameters.Nt); i++) {

      double T;

      switch(cooling_option){
        case 0:
          T = pow(10,T_array[i]);
          break;
        case 1:
          T = T_array[i];
          break;
      }

      mc_simulate(simulation_model,T);

      if(parameters.checkpoint_option){
        std::string save_loc {parameters.checkpoint_address + "structure_"
                              + std::to_string(i) + ".dat"};
        simulation_model.save_model_state(save_loc);
      }
      std::cout << "Energy at T = " << T << ": ";
      simulation_model.print_model_energy();
      std::cout << '\n' ;
    }
    std::string final_state_save_loc{parameters.final_structure_address +
                                     "final_structure.dat"};
    simulation_model.save_model_state(final_state_save_loc);
  }

  void mc::mc_simulate(model_space::model &simulation_model, double T){

    // Equilibrate the system for mcs_eq steps
    for (int step = 0; step < parameters.mcs_eq; step++) {
      simulation_model.update_model_system(T);
      simulation_model.update_model_records();
    }
    simulation_model.save_model_records(T);

    // Depending on the options in the mc_params structure, initialize the
    // containers that will store the MC averages
    simulation_model.initialize_model_averages();

    // Collect the averages
    for(int step=0;step<parameters.mcs_av;step++){
      simulation_model.update_model_system(T);
      simulation_model.update_model_averages(T);
    }

    /*Save averages to the files*/
    simulation_model.save_model_averages(T,parameters.mcs_av);
  }
}
