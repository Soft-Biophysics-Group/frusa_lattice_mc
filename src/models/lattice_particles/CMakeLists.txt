### Collect the header and source files to be compiled
set(HEADER_LATTICEPARTICLES
    ${INCLUDE_FRUSA_MODELS}/lattice_particles/lattice_particles_parameters.h
    ${INCLUDE_FRUSA_MODELS}/lattice_particles/lattice_particles_state.h
    ${INCLUDE_FRUSA_MODELS}/lattice_particles/lattice_particles_interactions.h
    ${INCLUDE_FRUSA_MODELS}/lattice_particles/lattice_particles_update.h
    ${INCLUDE_FRUSA_MODELS}/lattice_particles/lattice_particles_averages.h
    ${INCLUDE_FRUSA_THIRDPARTY}/json.hpp)

set(SOURCE_LATTICEPARTICLES
    lattice_particles_parameters.cc
    lattice_particles_state.cc
    lattice_particles_interactions.cc
    lattice_particles_update.cc
    lattice_particles_averages.cc
    )

### Create the lattice_particles library and include the header directories
add_library(lattice_particles_library
            ${SOURCE_LATTICEPARTICLES}
            ${HEADER_LATTICEPARTICLES})

target_include_directories(lattice_particles_library PUBLIC
                           ${INCLUDE_FRUSA_MODELS}/lattice_particles
                           ${INCLUDE_FRUSA_UTILITY}
                           ${INCLUDE_FRUSA_THIRDPARTY})

### Link the lattice_particles library to utils to make sure we have access to all of the
### helper functions
target_link_libraries(lattice_particles_library PRIVATE utils_library)

### Make sure that the compiler flags are the same for all libraries by
### linking them to compiler_flags
target_link_libraries(lattice_particles_library PRIVATE compiler_flags)
