add_library(vector_utils 
            vector_utils.cc
            ${HEADER_FRUSA_UTILITY})

target_include_directories(vector_utils PUBLIC 
                           ${PROJECT_SOURCE_DIR}/include/utility)

target_link_libraries(vector_utils PUBLIC compiler_flags)
