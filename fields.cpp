#include "fields.h" 

namespace model_space{
  /*
   * Definitions for the fields class
   */
  
  fields::fields(const model_data &model_data_1d) :
    N(model_data_1d.N),
    Np(model_data_1d.Np),
    k11(model_data_1d.k11),
    k12(model_data_1d.k12),
    k21(model_data_1d.k21),
    rng(model_data_1d.rng)
    {
    /*
     * Initialize the system of fields and calculate the initial energy
     */

    /*Calculate the total density*/
    psi_bar = (1.0*Np)/N;

    
    /*Initialize pseudorandom number distributions*/
    uniform_dist  = real_dist(0,1);
    particle_dist = int_dist(0,N-1);
    binary_dist   = int_dist(0,1);

    initialize();
    coupling_matrix = {{{k11,k21},{k12,k11}},
                       {{k11,k12},{k21,k11}}};

    energy = get_energy();

  }

  /*
   * Required public routines
   */

  void fields::print_state(){
    /*
     * Print the current state of the system
     */
    std::cout << "The current state of the system:\n\n";
    for(int i=0;i<N;i++){
      std::cout << concentrations[i][0] << " " << concentrations[i][1] << "\n";
    }
    std::cout << "\n";
  }

  void fields::save_state(std::string file_name, std::string address){
    /*
     * Save the current state of the system to a file
     */
    std::ofstream state_f;
    state_f.open(address+file_name);
    if(!state_f){
      std::cerr << "Could not open "+address+file_name << std::endl;
      exit(1);
    }

    for(int i=0;i<N;i++){
      state_f << std::setprecision(8) << concentrations[i][0] << " ";
      state_f << std::setprecision(8) << concentrations[i][1] << "\n";
    }
    state_f.close();
  }
  
  void fields::print_energy(){
    /*
     * Print the current energy of the system
     */
    std::cout << "energy = " << energy << "\n\n";
  }

  void fields::update_state(double T){
    /*
     * Update the state of the system using Metropolis algorithm
     */
    for(int i=0;i<N;i++){
      int site_index = particle_dist(rng);
      int update_type = binary_dist(rng);
      
      if(update_type==0){
        shift_local_density(site_index,T);
      }
      else{
        update_local_concentrations(site_index,T);
      }
    } 
  }
  
  void fields::initialize_averages(){
  }

  void fields::update_averages(double T){
  }

  void fields::save_averages(){
  }
  /*
   * Class-specific routines
   */

  void fields::initialize(){
    /*
     * Initialize a system of concentration fields on the lattice with 
     * uniform local density and random fractional concentrations 
     */
    
    /*Define and populate the concentration vector*/

    for(int i=0;i<N;i++){
      
      double psi_1 = psi_bar/2;//uniform_dist(rng)*psi_bar;
      double psi_2 = psi_bar-psi_1;

      concentrations.push_back({psi_1,psi_2});
    }
  }

  vec2i fields::get_neighbours(int r){
    
    /*
     * Extract the positions of the neighbours of a specified particle at 
     * position r
     */

    vec2i neighbours;

    int rm = (N+((r-1)%N))%N;
    int rp = (r+1)%N;

    neighbours.push_back({0,rm});
    neighbours.push_back({1,rp});

    return neighbours;
  }

  double fields::get_energy(){
    /*
     * Calculate total energy of the system
     */

    double en = 0;

    for(int i=0;i<N;i++){
      
      vec1d c_i = {concentrations[i][0],concentrations[i][1]};
      vec2i n = get_neighbours(i);

      for(int j=0;j<2;j++){
        int k = n[j][0];
        vec1d c_j = {concentrations[n[j][1]][0],concentrations[n[j][1]][1]};
 
        for(int a=0;a<2;a++){
          for(int b=0;b<2;b++){
            en+= coupling_matrix[k][a][b]*c_i[a]*c_j[b];
          }
        }
      }
    }
    return en/2;
  }

  void fields::shift_local_density(int site_index, double T){
    /*
     * Attempt to shift the local density from the selected lattice site
     * (site_index) to a new site with Metropolis acceptance rate
     */
    
    int donor_position = site_index;
    double donor_bound, acceptor_bound, bound_total, drho;
    int donor_type, acceptor_type;

    get_donor_bound(donor_position,donor_bound,donor_type);

    while(donor_bound==0.0){
      donor_position = (donor_position+1)%N;
      get_donor_bound(donor_position,donor_bound,donor_type);
    }

    /*Choose the lattice site where the density is transferred*/
    int acceptor_position = particle_dist(rng);

    if(acceptor_position==donor_position){
      acceptor_position = (acceptor_position+1)%N;
    }

    get_acceptor_bound(acceptor_position,acceptor_bound,acceptor_type);

    while(acceptor_bound==0.0){
      acceptor_position = (acceptor_position+1)%N;
      if(acceptor_position==donor_position){
        acceptor_position = (acceptor_position+1)%N;
      }
      get_acceptor_bound(acceptor_position,acceptor_bound,acceptor_type);
    }

    /*Total bound on local density transfer*/
    bound_total = std::min(donor_bound,acceptor_bound);

    /*Amount of density to transfer*/
    if(bound_total<1e-5){
      /*No need to shift tiny amounts of density (below threshold 1e-5), 
       *simply transfer all density from the site*/
      drho = bound_total;
    }
    else{
      drho = uniform_dist(rng)*bound_total;
    }

    /*Find the neighbours of the old and new positions*/
    vec2i n_d = get_neighbours(donor_position);
    vec2i n_a = get_neighbours(acceptor_position);

    /*Check if the old and new positions are first neighbours*/
    bool direct_neighbours;
    int dr;

    int rp1 = (donor_position+1)%N;
    int rm1 = (N+((donor_position-1)%N))%N;

    if(acceptor_position==rp1 or acceptor_position==rm1){
      direct_neighbours = true;
      dr = 1;
    }
    else{
      direct_neighbours = false;
      dr = 0;
    }

    /*Calculate the energy difference resulting from the position change*/
    double dE = 0;
    
    for(int j=0;j<2;j++){
      vec1d c_j_acceptor = {concentrations[n_a[j][1]][0],\
                            concentrations[n_a[j][1]][1]};
      
      vec1d c_j_donor = {concentrations[n_d[j][1]][0],\
                         concentrations[n_d[j][1]][1]};
 
      for(int k=0;k<2;k++){
        dE+= coupling_matrix[j][acceptor_type][k]*drho*c_j_acceptor[k];
        dE-= coupling_matrix[j][donor_type][k]*drho*c_j_donor[k];
      }
    }
    if(direct_neighbours){
      /*Adjust the energy difference in case the donor and acceptor sites 
        are nearest neighbours*/
      dE-= coupling_matrix[dr][donor_type][acceptor_type]*drho*drho;
    }

    /*Accept the new position using Metropolis rule*/
    if(dE<=0){
      concentrations[donor_position][donor_type]-= drho;
      concentrations[acceptor_position][acceptor_type]+= drho;
      energy+=dE;
    }
    else if(exp(-dE/T)>uniform_dist(rng)){
      concentrations[donor_position][donor_type]-= drho;
      concentrations[acceptor_position][acceptor_type]+= drho;
      energy+=dE;
    }
  }

  void fields::get_donor_bound(int donor_position, double &donor_bound,\
      int &donor_type){
    /*
     * Calculate the amount of local density that the donor site can transfer
     */

    vec1d donor_concentration = {concentrations[donor_position][0],\
                                 concentrations[donor_position][1]};

    /*Check if there is a single species type on the donor site*/
    if(donor_concentration[0]==0.0){
      donor_type = 1; 
    }
    else if(donor_concentration[1]==0.0){
      donor_type = 0;
    }
    else{
      donor_type = binary_dist(rng);
    }
    donor_bound = donor_concentration[donor_type];
  }

  void fields::get_acceptor_bound(int acceptor_position,\
      double &acceptor_bound, int &acceptor_type){
    /*
     * Calculate the amount of local density that the acceptor site can take
     */

    vec1d acceptor_concentration = {concentrations[acceptor_position][0],\
                                    concentrations[acceptor_position][1]};

    /*Check if the site is full*/
    double acceptor_local_density = acceptor_concentration[0]+\
                                    acceptor_concentration[1];
    if(1.0-acceptor_local_density<1e-5){
      acceptor_bound = 0.0; 
    }
    else{
      acceptor_type = binary_dist(rng); 
      acceptor_bound = std::min(1-acceptor_concentration[acceptor_type],\
                                1-acceptor_local_density);
    }
  }

  void fields::update_local_concentrations(int site_index, double T){
    /* 
     * Attempt to re-distribute the fractional concentrations on the selected
     * lattice site (site_index) with Metropolis acceptance rate
     */
    
    /*Determine which fractional concentration is reduced and which one is 
      increased*/
    int donor_position=site_index;
    double donor_bound, drho;
    int donor_type, acceptor_type;

    get_donor_bound(donor_position,donor_bound,donor_type);
    
    while(donor_bound==0.0){
      donor_position = (donor_position+1)%N;
      get_donor_bound(donor_position,donor_bound,donor_type);
    }

    acceptor_type = 1-donor_type;

    /*Amount of concentration to convert*/
    if(donor_bound<1e-5){
      /*No need to shift tiny amounts of density (below threshold 1e-5), 
       *simply convert all concentration to the other type*/
      drho = donor_bound;
    }
    else{
      drho = uniform_dist(rng)*donor_bound;
    }

    /*Find the neighbours of the selected site*/
    vec2i n = get_neighbours(donor_position);

    /*Calculate the energy difference resulting from the position change*/
    double dE = 0;
    
    for(int j=0;j<2;j++){
      vec1d c_j = {concentrations[n[j][1]][0],\
                   concentrations[n[j][1]][1]};
 
      for(int k=0;k<2;k++){
        dE+= coupling_matrix[j][acceptor_type][k]*drho*c_j[k];
        dE-= coupling_matrix[j][donor_type][k]*drho*c_j[k];
      }
    }

    /*Accept the new position using Metropolis rule*/
    if(dE<=0){
      concentrations[donor_position][donor_type]-= drho;
      concentrations[donor_position][acceptor_type]+= drho;
      energy+=dE;
    }
    else if(exp(-dE/T)>uniform_dist(rng)){
      concentrations[donor_position][donor_type]-= drho;
      concentrations[donor_position][acceptor_type]+= drho;
      energy+=dE;
    }
  }

  void fields::update_psi(vec2d &psi){
    /*
     * Update the density vector
     */

    /*for(int i=0;i<Np;i++){
      psi[positions[i]][orientations[i]-1] += 1;
    }*/
  }
}

