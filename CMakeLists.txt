# To find the cmake version use 
# $ cmake --version
cmake_minimum_required(VERSION 3.16...3.22)

# Set required compiler to 2011 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the name of the project
project(fusa_mc VERSION 1.0.0
  DESCRIPTION 
  "Monte Carlo simulation package for models of frustrated self-assembly"
  LANGUAGES CXX)

# Add the model library
add_subdirectory(src/models)

# Set directory with the input files
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input/model_params.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input/mc_params.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Add an executable with the above sources
add_executable(frusa_mc app/main.cpp)

# Link the model library to the executable
target_link_libraries(frusa_mc PUBLIC model)

# Set the directory with the header files
target_include_directories(frusa_mc PRIVATE 
  ${PROJECT_SOURCE_DIR}/include
)

