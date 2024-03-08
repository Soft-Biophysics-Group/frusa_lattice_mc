add_library(vector_utils 
            vector_utils.cc
            ${HEADER_FRUSA_UTILITY})

target_include_directories(vector_utils PUBLIC 
                           ${INCLUDE_FRUSA_UTILITY})

target_link_libraries(vector_utils PUBLIC compiler_flags)
