set(SOURCE_FIELDS 
    fields_state.cc
    fields_interactions.cc
    fields_chain.cc)

add_library(fields_library ${SOURCE_FIELDS})

target_include_directories(fields_library PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/models/fields
                           ${INCLUDE_FRUSA_UTILITY}
                           ${INCLUDE_FRUSA_THIRDPARTY})

target_link_libraries(fields_library PUBLIC compiler_flags)
target_link_libraries(fields_library PRIVATE vector_utils)
