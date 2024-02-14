#include "mc_routines.h"

namespace simulation{
  template <class model>
  mc<model>::mc(model m, const mc_data &mc_data_sim) :
    lattice_system(m),
    mcs_eq(mc_data_sim.mcs_eq),
    mcs_av(mc_data_sim.mcs_av),
    cooling_schedule(mc_data_sim.cooling_schedule),
    Ti(mc_data_sim.Ti),
    Tf(mc_data_sim.Tf),
    Nt(mc_data_sim.Nt),
    checkpoint(mc_data_sim.checkpoint)
    {
    /*
     * Initialize the Monte Carlo simulation on a selected model
     */

    /*Define the array of temperatures for the annealing*/
    dT = (Tf-Ti)/Nt;
        
    for(int i=0;i<Nt;i++){
      T_array.push_back(Ti+i*dT);
    }

    /*If the checkpoint option is selected, define the corresponding output
      address*/
    if(checkpoint){
      checkpoint_address = mc_data_sim.checkpoint_address;
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

      switch(cooling_schedule){
        /*Depending on the cooling schedule, define the temperature for 
         * Metropolis probabilities*/
        case "exponential":
          T = pow(10,T_array[i]);
          break;
        case "linear" :
          T = T_array[i];
          break;
      }
      mc_simulate(T);
      if(checkpoint){
        lattice_system.save_state(checkpoint_address+"structure_"+\
                                  std::to_string(i)+".dat");
      }
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

    /*Depending on the options in the mc_data structure, initialize the 
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
