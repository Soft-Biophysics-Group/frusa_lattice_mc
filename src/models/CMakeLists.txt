### Create the main model library

add_library(model_library model.cc
            ${HEADER_FRUSA_MODELS})

target_include_directories(model_library PUBLIC
                           ${INCLUDE_FRUSA_MODELS}
                           ${INCLUDE_FRUSA_THIRDPARTY})

target_link_libraries(model_library PRIVATE compiler_flags)

### Create a compile option to select user-defined model library

if("${MODEL_TYPE}" STREQUAL "default")
  target_compile_definitions(model_library PUBLIC "DEFAULT")

  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/default)

  target_include_directories(model_library PUBLIC
                             ${INCLUDE_FRUSA_MODELS}/default)

  target_link_libraries(model_library PRIVATE default_library)


elseif("${MODEL_TYPE}" STREQUAL "fields_chain" OR
   "${MODEL_TYPE}" STREQUAL "fields_square" OR
   "${MODEL_TYPE}" STREQUAL "fields_hexagonal")

  target_compile_definitions(model_library PUBLIC "FIELDS")

  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fields)

  target_include_directories(model_library PUBLIC
                             ${INCLUDE_FRUSA_MODELS}/fields)

  target_link_libraries(model_library PRIVATE fields_library)

elseif("${MODEL_TYPE}" STREQUAL "lattice_hexagonal")
  target_compile_definitions(model_library PUBLIC "LATTICEPARTICLES")

  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lattice_particles)

  target_include_directories(model_library PUBLIC
                             ${INCLUDE_FRUSA_MODELS}/lattice_particles)

  target_link_libraries(model_library PRIVATE lattice_particles_library)
endif()
