// Copyright (c) 2024 Soft Biophysics Group LPTMS
// Part of frusa_mc, released under BSD 3-Clause License.

#include "io_utils.h"

namespace io_space{

  /*
   * Routines for reading data from files
   */

  void read_vector(vec1i &array, int N, std::string location){

    std::ifstream array_f;
    array_f.open(location);

    if(!array_f){
      std::cerr << "Unable to open file "+location;
      exit(1);
    }

    for(int i=0;i<N;i++){

      int a_i;
      array_f >> a_i;

      array.push_back(a_i);
    }
    array_f.close();
  }

  void read_vector(vec1d &array, int N, std::string location){

    std::ifstream array_f;
    array_f.open(location);

    if(!array_f){
      std::cerr << "Unable to open file "+location;
      exit(1);
    }

    for(int i=0;i<N;i++){

      double a_i;
      array_f >> a_i;

      array.push_back(a_i);
    }
    array_f.close();
  }

  void read_vector(vec2i &array, int N1, int N2, std::string location){

    std::ifstream array_f;
    array_f.open(location);

    if(!array_f){
      std::cerr << "Unable to open file "+location;
      exit(1);
    }

    for(int i=0;i<N1;i++){

      vec1i a_i;

      for(int j=0;j<N2;j++){

        int a_ij;
        array_f >> a_ij;

        a_i.push_back(a_ij);
      }
      array.push_back(a_i);
    }
    array_f.close();
  }

  void read_vector(vec2d &array, int N1, int N2, std::string location){

    std::ifstream array_f;
    array_f.open(location);

    if(!array_f){
      std::cerr << "Unable to open file "+location;
      exit(1);
    }

    for(int i=0;i<N1;i++){

      vec1d a_i;

      for(int j=0;j<N2;j++){

        double a_ij;
        array_f >> a_ij;

        a_i.push_back(a_ij);
      }
      array.push_back(a_i);
    }
    array_f.close();
  }

  /*
   * Routines for writing data to files
   */

  void save_vector(vec1i &array, int N, std::string location){

    std::ofstream array_f;
    array_f.open(location);

    if(!array_f){
      std::cerr << "Unable to open file "+location;
      exit(1);
    }

    for (std::size_t i = 0; i < static_cast<std::size_t>(N); i++) {
      array_f << array[i] << "\n";
    }
    array_f.close();
  }

  void save_vector(vec1d &array, int N, std::string location){

    std::ofstream array_f;
    array_f.open(location);

    if(!array_f){
      std::cerr << "Unable to open file "+location;
      exit(1);
    }

    for(std::size_t i=0;i<static_cast<std::size_t>(N);i++){
      array_f << std::setprecision(8) << array[i] << "\n";
    }
    array_f.close();
  }

  void save_vector(vec2i &array, int N1, int N2, std::string location){

    std::ofstream array_f;
    array_f.open(location);

    if(!array_f){
      std::cerr << "Unable to open file "+location;
      exit(1);
    }

    for(std::size_t i=0;i<static_cast<std::size_t>(N1);i++){
      for(std::size_t j=0;j<static_cast<std::size_t>(N2);j++){
        array_f << array[i][j] << " ";
      }
      array_f << "\n";
    }
    array_f.close();
  }

  void save_vector(vec2d &array, int N1, int N2, std::string location){

    std::ofstream array_f;
    array_f.open(location);

    if(!array_f){
      std::cerr << "Unable to open file "+location;
      exit(1);
    }

    for(std::size_t i=0;i<static_cast<std::size_t>(N1);i++){
      for(std::size_t j=0;j<static_cast<std::size_t>(N2);j++){
        array_f << std::setprecision(8) << array[i][j] << " ";
      }
      array_f << "\n";
    }
    array_f.close();
  }
}
