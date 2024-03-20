#include "fields_hexagonal.h"

namespace fields_space{
  
  vec3d get_coupling_matrix_hexagonal(vec1d couplings){
    
    vec3d coupling_matrix;
    
    // Read the coupling matrix for the bond along the +x axis from the input
    // file. The structure of couplings is assumed to be a flattened version of
    // the coupling_matrix.
    vec2d array_0;
    for(int a=0;a<6;a++){
      vec1d row;
      for(int b=0;b<6;b++){
        row.push_back(couplings[1+a+6*b]);
      }
      array_0.push_back(row);
    }

    coupling_matrix.push_back(array_0);

    for(int k=1;k<6;k++){
      vec2d array_k;
      for(int a=0;a<6;a++){
        vec1d row;
        for(int b=0;b<6;b++){
          int ap = array_space::mod(a-k,6);
          int bp = array_space::mod(a+k,6);
          row.push_back(coupling_matrix[0][ap][bp]);
        }
        array_k.push_back(row);
      }
      coupling_matrix.push_back(array_k);
    }
    
    return coupling_matrix; 
  }

  vec1i get_neighbours_hexagonal(int r, int Lx, int Ly){
    
    int i,j;

    array_space::r_to_ij(r,i,j,Lx,Ly);

    int im = array_space::mod(i-1,Lx);
    int ip = array_space::mod(i+1,Lx);
    int jm = array_space::mod(j-1,Ly);
    int jp = array_space::mod(j+1,Ly);

    int r0,r1,r2,r3;
    array_space::ij_to_r(r0,ip,j ,Lx,Ly);
    array_space::ij_to_r(r1,ip,jp,Lx,Ly);
    array_space::ij_to_r(r2,i ,jp,Lx,Ly);
    array_space::ij_to_r(r3,im,j ,Lx,Ly);
    array_space::ij_to_r(r4,im,jm,Lx,Ly);
    array_space::ij_to_r(r5,i ,jm,Lx,Ly);

    vec1i neighbours = {r0,r1,r2,r3,r4,r5};
    
    return neighbours;
  }
  
  int get_bond_direction_hexagonal(int r1, int r2, int Lx, int Ly){

    int i1,j1;
    array_space::r_to_ij(r1,i1,j1,Lx,Ly);

    int dr;

    int i1m = array_space::mod(i1-1,Lx);
    int i1p = array_space::mod(i1+1,Lx);
    int j1m = array_space::mod(j1-1,Ly);
    int j1p = array_space::mod(j1+1,Ly);

    int r1_0,r1_1,r1_2,r1_3;
    array_space::ij_to_r(r1_0,i1p,j1 ,Lx,Ly);
    array_space::ij_to_r(r1_1,i1p,j1p,Lx,Ly);
    array_space::ij_to_r(r1_2,i1 ,j1p,Lx,Ly);
    array_space::ij_to_r(r1_3,i1m,j1 ,Lx,Ly);
    array_space::ij_to_r(r1_4,i1m,j1m,Lx,Ly);
    array_space::ij_to_r(r1_5,i1 ,j1m,Lx,Ly);

    if(r2==r1_0){
      dr = 0;
    }
    else if(r2==r1_1){
      dr = 1;
    }
    else if(r2==r1_2){
      dr = 2;
    }
    else if(r2==r1_3){
      dr = 3;
    }
    else if(r2==r1_4){
      dr = 4;
    }
    else if(r2==r1_5){
      dr = 5;
    }
    return dr;
  }
}
