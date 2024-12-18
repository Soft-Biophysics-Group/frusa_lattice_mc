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
add_compile_options(-Weffc++ -Wextra -Wconversion -Wsign-conversion -Werror)
set (CMAKE_CXX_STANDARD 20)

################################################################################


### Gather headers and sources #################################################

### Header Groups
set(HEADER_FRUSA_ENGINE)
set(HEADER_FRUSA_MODELS)
set(HEADER_FRUSA_GEOMETRY)
set(HEADER_FRUSA_UTILITY)
set(HEADER_FRUSA_THIRDPARTY)

### Header directories
set(INCLUDE_FRUSA_ENGINE ${PROJECT_SOURCE_DIR}/include/engine)
set(INCLUDE_FRUSA_MODELS ${PROJECT_SOURCE_DIR}/include/models)
set(INCLUDE_FRUSA_GEOMETRY ${PROJECT_SOURCE_DIR}/../include/geometry)
set(INCLUDE_FRUSA_UTILITY ${PROJECT_SOURCE_DIR}/../include/utility)
set(INCLUDE_FRUSA_THIRDPARTY ${PROJECT_SOURCE_DIR}/../include/thirdparty)

### Propagate variables to project subdirectories
add_subdirectory(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/engine)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/models)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/utility)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/geometry)

###############################################################################

### Create the executable

add_executable(frusa_lattice_mc_test test_geometry.cc
               ${HEADER_FRUSA_GEOMETRY})

### Link the libraries to the executable
target_link_libraries(frusa_lattice_mc_test PUBLIC compiler_flags)
target_link_libraries(frusa_lattice_mc_test PUBLIC geometry_library)
target_link_libraries(frusa_lattice_mc_test PUBLIC models_library)
target_link_libraries(frusa_lattice_mc_test PUBLIC particles_library)
target_link_libraries(frusa_lattice_mc_test PUBLIC mc_library)
target_link_libraries(frusa_lattice_mc_test PUBLIC utils_library)
