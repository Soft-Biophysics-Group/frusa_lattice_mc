set(HEADER_FRUSA_MODELS
    ${HEADER_FRUSA_MODELS}
    ${CMAKE_CURRENT_SOURCE_DIR}/lattice_particles_state.h
    ${CMAKE_CURRENT_SOURCE_DIR}/lattice_particles_interactions.h
    PARENT_SCOPE)

if("${MODEL_TYPE}" STREQUAL "lattice_hexagonal")

  set(HEADER_FRUSA_MODELS
      ${HEADER_FRUSA_MODELS}
      ${CMAKE_CURRENT_SOURCE_DIR}/lattice_particles_hexagonal.h
      PARENT_SCOPE)

endif()
