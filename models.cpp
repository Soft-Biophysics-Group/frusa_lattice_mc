#include "models.h" 

namespace simulation{
  /*
   * Definitions for the particle class
   */
  
  particles::particles(const model_data &model_data_1d){
    /*
     * Initialize the system of particles and calculate the initial energy
     */

    N   = model_data_1d.N;
    Np  = model_data_1d.Np;
    k11 = model_data_1d.k11;
    k12 = model_data_1d.k12;
    k21 = model_data_1d.k21;

    initialize();
    coupling_matrix = {{{k11,k21},{k12,k11}},
                       {{k11,k12},{k21,k11}}};

    energy = get_energy();

  }

  /*
   * Required public routines
   */

  void particles::print_state(){
    /*
     * Print the current state of the system
     */
    for(int i=0;i<Np;i++){
      std::cout << state[i][0] << ", " << state[i][1] << "\n";
    }
  }

  void particles::print_energy(){
    /*
     * Print the current energy of the system
     */
    std::cout << energy << "\n";
  }

  /*
   * Class-specific routines
   */

  void particles::initialize(){
    /*
     * Initialize a system of particles at random locations and with
     * random orientations 
     */
    EngineType rng(dev());
    
    /*Define a re-shuffled vector with all positions on a 1D lattice*/
    real_dist u(0,1);   

    vec1i all_positions(N);

    for(int i=0;i<N;i++){
      
      all_positions[i] = i;
    
    }

    std::shuffle(std::begin(all_positions), std::end(all_positions), rng);

    /*Define and populate the state vector*/

    for(int i=0;i<Np;i++){
       
      if(u(rng)>0.5){
        state.push_back({all_positions[i],1});
      }
      else{
        state.push_back({all_positions[i],2});
      }
    }
  }

  vec1i particles::get_neighbours(int r){
    
    /*
     * Extract the indices (location in the state vector) of the neighbours
     * of a specified particle at position r
     */

    vec1i neighbours={-1,-1};

    int rm = (N+((r-1)%N))%N;
    int rp = (r+1)%N;

    for(int i=0;i<Np;i++){
      if(state[i][0]==rm){
        neighbours[0] = i;
      }
      else if(state[i][0]==rp){
        neighbours[1] = i;
      }
    }

    return neighbours;
  }

  double particles::get_energy(){
    /*
     * Calculate total energy of the system
     */

    double en = 0;

    for(int i=0;i<Np;i++){
      
      vec1i n = get_neighbours(state[i][0]);

      for(int j=0;j<2;j++){
      
        if(n[j]!=-1){
          en+= coupling_matrix[j][ state[i][1]-1][ state[ n[j] ][1]-1];
        }
      }
    }
    return en/2;
  }

  void particles::update_psi(vec2d &psi){
    /*
     * Update the density vector
     */

    for(int i=0;i<Np;i++){
      psi[state[i][0]][state[i][1]-1] += 1;
    }
  }
}

