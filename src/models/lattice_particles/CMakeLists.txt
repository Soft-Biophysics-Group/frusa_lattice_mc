### Collect the header and source files to be compiled
set(HEADER_LATTICEPARTICLES
    ${INCLUDE_FRUSA_MODELS}/lattice_particles/lattice_particles_parameters.h
    ${INCLUDE_FRUSA_MODELS}/lattice_particles/lattice_particles_state.h
    #${INCLUDE_FRUSA_MODELS}/lattice_particles/lattice_particles_interactions.h
    #${INCLUDE_FRUSA_MODELS}/lattice_particles/lattice_particles_update.h
    #${INCLUDE_FRUSA_MODELS}/lattice_particles/lattice_particles_averages.h
    ${INCLUDE_FRUSA_THIRDPARTY}/json.hpp)

set(SOURCE_LATTICEPARTICLES
    lattice_particles_parameters.cc
    lattice_particles_state.cc
    #lattice_particles_interactions.cc
    #lattice_particles_update.cc
    #lattice_particles_averages.cc
    )

### Create the lattice_particles library and include the header directories
add_library(lattice_particles_library
            ${SOURCE_LATTICEPARTICLES}
            ${HEADER_LATTICEPARTICLES})

if("${MODEL_TYPE}" STREQUAL "lattice_hexagonal")

    target_compile_definitions(fields_library PRIVATE "LATTICEPARTICLES_HEXAGONAL")

    add_library(user_library STATIC lattice_particles_hexagonal.cc)

endif()

#target_include_directories(user_library PUBLIC
                           #${INCLUDE_FRUSA_MODELS}/lattice_particles
                           #${INCLUDE_FRUSA_UTILITY})

#target_link_libraries(user_library PRIVATE compiler_flags)
#target_link_libraries(user_library PRIVATE vector_utils)

#target_link_libraries(lattice_particles_library PRIVATE user_library)

target_include_directories(lattice_particles_library PUBLIC
                           ${INCLUDE_FRUSA_MODELS}/lattice_particles
                           ${INCLUDE_FRUSA_UTILITY}
                           ${INCLUDE_FRUSA_THIRDPARTY})

target_link_libraries(lattice_particles_library PRIVATE compiler_flags)
target_link_libraries(lattice_particles_library PRIVATE vector_utils)
