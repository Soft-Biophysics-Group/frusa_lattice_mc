set(HEADER_DEFAULT
    ${INCLUDE_FRUSA_MODELS}/default/default_parameters.h
    ${INCLUDE_FRUSA_MODELS}/default/default_state.h
    ${INCLUDE_FRUSA_MODELS}/default/default_interactions.h
    ${INCLUDE_FRUSA_MODELS}/default/default_update.h
    ${INCLUDE_FRUSA_THIRDPARTY}/json.hpp)

set(SOURCE_DEFAULT
    default_parameters.cc
    default_state.cc
    default_interactions.cc
    default_update.cc)

add_library(default_library 
            ${SOURCE_DEFAULT}
            ${HEADER_DEFAULT})

target_include_directories(default_library PUBLIC
                           ${INCLUDE_FRUSA_MODELS}/default
                           ${INCLUDE_FRUSA_UTILITY}
                           ${INCLUDE_FRUSA_THIRDPARTY})

target_link_libraries(default_library PRIVATE compiler_flags)
target_link_libraries(default_library PRIVATE vector_utils)

