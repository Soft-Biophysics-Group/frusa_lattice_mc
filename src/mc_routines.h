#ifndef MC_HEADER_H
#define MC_HEADER_H

#include "utils.h"

namespace simulation_space{
  
  /*Options for the cooling schedule*/
  enum cooling_option {exponential, linear};

  /*Monte Carlo parameters*/
  struct mc_data{
    int mcs_eq;
    int mcs_av;
    double Ti;
    double Tf;
    int Nt;
    cooling_option cooling_schedule;
    bool checkpoint;
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
      cooling_option cooling_schedule;

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
      bool checkpoint;

      /*String to store the checkpoint output address*/
      std::string checkpoint_address;

    public:

      /*Class constructor*/
      mc(model, const struct mc_data &);

      /*MC annealing*/
      void t_scan();

      /*MC simmulation at a fixed temperature*/
      void mc_simulate(double);
 };
  
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
        case exponential:
          T = pow(10,T_array[i]);
          break;
        case linear:
          T = T_array[i];
          break;
      }
      mc_simulate(T);
      if(checkpoint){
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

#endif
