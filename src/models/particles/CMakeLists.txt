### Collect the header and source files to be compiled
set(HEADER_PARTICLES
    ${INCLUDE_FRUSA_MODELS}/particles/particles_parameters.h
    ${INCLUDE_FRUSA_MODELS}/particles/particles_state.h
    ${INCLUDE_FRUSA_MODELS}/particles/particles_interactions.h
    ${INCLUDE_FRUSA_MODELS}/particles/particles_update.h
    ${INCLUDE_FRUSA_MODELS}/particles/particles_averages.h
    ${INCLUDE_FRUSA_THIRDPARTY}/json.hpp)

set(SOURCE_PARTICLES
    particles_parameters.cc
    particles_state.cc
    particles_interactions.cc
    particles_update.cc
    particles_averages.cc
    )

### Create the particles library and include the header directories
add_library(particles_library
            ${SOURCE_PARTICLES}
            ${HEADER_PARTICLES})

if("${MODEL_TYPE}" STREQUAL "particles_hexagonal")

    target_compile_definitions(particles_library PRIVATE "PARTICLES_HEXAGONAL")

    add_library(user_library STATIC particles_hexagonal.cc)

endif()

target_include_directories(user_library PUBLIC
                           ${INCLUDE_FRUSA_MODELS}/particles
                           ${INCLUDE_FRUSA_GEOMETRY}
                           ${INCLUDE_FRUSA_UTILITY}
                           ${INCLUDE_FRUSA_THIRDPARTY})

target_link_libraries(user_library PRIVATE compiler_flags)
target_link_libraries(user_library PRIVATE utils_library)
target_link_libraries(user_library PRIVATE geometry_library)

target_link_libraries(particles_library PRIVATE user_library)

target_include_directories(particles_library PUBLIC
                           ${INCLUDE_FRUSA_MODELS}/particles
                           ${INCLUDE_FRUSA_UTILITY}
                           ${INCLUDE_FRUSA_GEOMETRY}
                           ${INCLUDE_FRUSA_THIRDPARTY})

target_link_libraries(particles_library PRIVATE compiler_flags)
target_link_libraries(particles_library PRIVATE utils_library)
target_link_libraries(particles_library PRIVATE geometry_library)
