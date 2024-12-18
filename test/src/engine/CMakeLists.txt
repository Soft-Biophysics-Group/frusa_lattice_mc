### Create the mc library

set(SOURCE_FRUSA_ENGINE
    ${CMAKE_CURRENT_SOURCE_DIR}/mc_routines.cc)

add_library(mc_library
            ${SOURCE_FRUSA_ENGINE}
            ${HEADER_FRUSA_ENGINE}
            ${HEADER_FRUSA_MODELS}
            ${HEADER_FRUSA_UTILITY}
            ${HEADER_FRUSA_THIRDPARTY})

target_include_directories(mc_library PUBLIC
                           ${INCLUDE_FRUSA_ENGINE}
                           ${INCLUDE_FRUSA_MODELS}
                           ${INCLUDE_FRUSA_GEOMETRY}
                           ${INCLUDE_FRUSA_UTILITY}
                           ${INCLUDE_FRUSA_THIRDPARTY})

target_link_libraries(mc_library PRIVATE compiler_flags)
target_link_libraries(mc_library PRIVATE model_library)
target_link_libraries(mc_library PRIVATE utils_library)

