#ifndef MC_HEADER_H
#define MC_HEADER_H

#include "utils.h"

namespace simulation_space{
  
  /*Monte Carlo parameters*/
  struct mc_params{
    mc_params();
    int mcs_eq;
    int mcs_av;
    double Ti;
    double Tf;
    int Nt;
    std::string cooling_schedule;
    bool checkpoint_option;
    std::string checkpoint_address;
  };
     
  template <class model> 
  class mc {
    private:
      /*
       * Simulation parameters
       */

      /*Model studied with MC*/

      model lattice_system;

      /*Number of MC steps for equilibration and averaging*/
      int mcs_eq, mcs_av;

      /*Type of cooling schedule*/
      std::string cooling_schedule;
      int cooling_option;

      /*Array of temperature values*/
      vec1d T_array;

      /*Initial/final temperatures, step size*/
      double Ti, Tf, dT;

      /*Number of temperature steps*/
      int Nt;

      /*
       * Options for data collections
       */

      /*Save state configurations at intermediate temperatures?*/
      bool checkpoint_option;

      /*String to store the checkpoint output address*/
      std::string checkpoint_address;

    public:

      /*Class constructor*/
      mc(model, const struct mc_params &);

      /*MC annealing*/
      void t_scan();

      /*MC simmulation at a fixed temperature*/
      void mc_simulate(double);
 };
 
  mc_params::mc_params(){
    /*
     * Populate the struct using the input JSON file
     */
    
    std::ifstream mc_f;
    mc_f.open("./mc_params.json");
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
  }

  template <class model>
  mc<model>::mc(model m, const mc_params &mc_params_sim) :
    lattice_system(m),
    mcs_eq(mc_params_sim.mcs_eq),
    mcs_av(mc_params_sim.mcs_av),
    cooling_schedule(mc_params_sim.cooling_schedule),
    Ti(mc_params_sim.Ti),
    Tf(mc_params_sim.Tf),
    Nt(mc_params_sim.Nt),
    checkpoint_option(mc_params_sim.checkpoint_option)
    {
    /*
     * Initialize the Monte Carlo simulation on a selected model
     */

    /*Map the cooling option on the integer variable*/
    try{
      if(cooling_schedule=="exponential"){
        cooling_option = 0;
      }
      else if(cooling_schedule=="linear"){
        cooling_option = 1;
      }
      else{
        throw cooling_schedule;
      }
    }

    catch(std::string cooling_schedule){
      std::cout << cooling_schedule << ": Incorrect cooling option!\n";
      exit(1);
    }

    /*Define the array of temperatures for the annealing*/
    dT = (Tf-Ti)/Nt;
        
    for(int i=0;i<Nt;i++){
      T_array.push_back(Ti+i*dT);
    }

    /*If the checkpoint option is selected, define the corresponding output
      address*/
    if(checkpoint_option){
      checkpoint_address = mc_params_sim.checkpoint_address;
    }

  }

  template <class model>
  void mc<model>::t_scan(){
    /*
     * Performs MC annealing according to the temperature range defined in
     * T_array
     */
   
    for(int i=0;i<Nt;i++){
      
      double T;

      switch(cooling_option){
        /*Depending on the cooling schedule, define the temperature for 
         * Metropolis probabilities*/
        case 0:
          T = pow(10,T_array[i]);
          break;
        case 1:
          T = T_array[i];
          break;
      }
      mc_simulate(T);
      if(checkpoint_option){
        lattice_system.save_state("structure_"+std::to_string(i)+".dat",\
                                  checkpoint_address);
      }
      std::cout << "Energy at T = " << T << ":\n";
      lattice_system.print_energy();
    }
  }

  template <class model>
  void mc<model>::mc_simulate(double T){
    /*
     * Performs MC simulation at a fixed temperature T
     */

    /*Equilibrate the system for mcs_eq steps*/
    for(int step=0;step<mcs_eq;step++){  
      lattice_system.update_state(T);
    }

    /*Depending on the options in the mc_params structure, initialize the 
     *containers that will store the MC averages*/
    lattice_system.initialize_averages();

    /*Collect the averages*/
    for(int step=0;step<mcs_av;step++){
      lattice_system.update_state(T);
      lattice_system.update_averages(T);
    }

    /*Save averages to the files*/
    lattice_system.save_averages();
  }
}

#endif
