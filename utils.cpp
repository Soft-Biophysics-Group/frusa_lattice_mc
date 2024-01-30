#include "utils.h" 

namespace lattice_system{
  particles::particles(int N, int Np, double k11, double k12, double k21){
    /*
     * Initialize the system of particles and calculate the initial energy
     */

    initialize(N,Np);
    coupling_matrix = {{{k11,k21},{k12,k11}},
                       {{k11,k12},{k21,k11}}};

    E = get_energy(N,Np);

  }

  int particles::initialize(int N, int Np){
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
    return 0;
  }

  vec1i particles::get_neighbours(int r, int N, int Np){
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

  double particles::get_energy(int N, int Np){
    /*
     * Calculate total energy of the system
     */

    double energy = 0;

    for(int i=0;i<Np;i++){
      
      vec1i n = get_neighbours(state[i][0],N,Np);

      for(int j=0;j<2;j++){
      
        if(n[j]!=-1){
          energy+= coupling_matrix[j][ state[i][1]-1][ state[ n[j] ][1]-1];
        }
      }
    }
    return energy/2;
  }

  int particles::update_psi(vec2d &psi, int Np){
    /*
     * Update the density vector
     */

    for(int i=0;i<Np;i++){
      psi[state[i][0]][state[i][1]-1] += 1;
    }
    
    return 0;
  }
}

