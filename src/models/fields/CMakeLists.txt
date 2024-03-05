add_library(fields_library fields_state.cc)

target_include_directories(fields_library PUBLIC
                           "${PROJECT_SOURCE_DIR}/include/models/fields"
                           "${PROJECT_SOURCE_DIR}/include/thirdparty")

target_link_libraries(fields_library PUBLIC compiler_flags)
target_link_libraries(fields_library PRIVATE vector_utils)
