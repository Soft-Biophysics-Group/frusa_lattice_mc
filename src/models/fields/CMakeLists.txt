set(SOURCE_FIELDS 
    ${CMAKE_CURRENT_SOURCE_DIR}/fields_state.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/fields_interactions.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/fields_chain.cc)

add_library(fields_library ${SOURCE_FIELDS})

target_include_directories(fields_library PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/models/fields
                           ${PROJECT_SOURCE_DIR}/include/utility
                           ${PROJECT_SOURCE_DIR}/include/thirdparty)

target_link_libraries(fields_library PUBLIC compiler_flags)
target_link_libraries(fields_library PRIVATE vector_utils)
