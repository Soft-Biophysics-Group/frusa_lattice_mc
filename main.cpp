#include "utils.h"

int main(){

  int N=20;
  int Np=5;

  double k11 = -1;
  double k12 = 0;
  double k21 = 0;

  lattice_system::particles test(N,Np,k11,k12,k21);

  for(int i=0;i<Np;i++){
    std::cout << test.state[i][0] << ", " << test.state[i][1] << "\n";
  }

  std::cout << "\nE = " << test.E << "\n";

  vec2d psi(N,vec1d(2));

  test.update_psi(psi,Np);

  for(int i=0;i<N;i++){
    std::cout << psi[i][0] << "  " << psi[i][1] << "\n";
  }

}
