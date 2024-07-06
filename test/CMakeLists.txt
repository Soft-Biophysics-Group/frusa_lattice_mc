### CMake Version ##############################################################

cmake_minimum_required(VERSION 3.16...3.22)

################################################################################


### Create frusa_mc project ####################################################

project(frusa_lattice_mc_test 
        LANGUAGES CXX)

################################################################################


### Compiler options ###########################################################

add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_11)

################################################################################


### Gather headers and sources #################################################

### Header Groups
set(HEADER_FRUSA_GEOMETRY)
set(HEADER_FRUSA_UTILITY)
set(HEADER_FRUSA_THIRDPARTY)

### Header directories
set(INCLUDE_FRUSA_GEOMETRY ${PROJECT_SOURCE_DIR}/../include/geometry)
set(INCLUDE_FRUSA_UTILITY ${PROJECT_SOURCE_DIR}/../include/utility)
set(INCLUDE_FRUSA_THIRDPARTY ${PROJECT_SOURCE_DIR}/../include/thirdparty)

### Propagate variables to project subdirectories
add_subdirectory(${PROJECT_SOURCE_DIR}/../include/geometry)
add_subdirectory(${PROJECT_SOURCE_DIR}/../include/utility)
add_subdirectory(${PROJECT_SOURCE_DIR}/../src/geometry)
add_subdirectory(${PROJECT_SOURCE_DIR}/../src/utility)

###############################################################################

### Create the executable

add_executable(frusa_lattice_mc_test main.cc
               ${HEADER_FRUSA_GEOMETRY})

### Link the libraries to the executable
target_link_libraries(frusa_lattice_mc_test PUBLIC compiler_flags)
target_link_libraries(frusa_lattice_mc_test PUBLIC geometry_library)
