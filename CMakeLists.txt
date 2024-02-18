# To find the cmake version use 
# $ cmake --version
cmake_minimum_required(VERSION 3.22.1)

# Set the name of the project
project(fusa_mc VERSION 1.0.0
  DESCRIPTION 
  "Monte Carlo simulation package for models of frustrated self-assembly"
  LANGUAGES CXX)

# Create a sources variable with a link to the cpp files to compile
set(SOURCES 
  src/fields.cpp
  src/particles.cpp
  src/main.cpp
)

# Add an executable with the above sources
add_executable(frusa_mc ${SOURCES})

# Set the directory with the header files
target_include_directories(frusa_mc PRIVATE 
  ${PROJECT_SOURCE_DIR}/include
)
