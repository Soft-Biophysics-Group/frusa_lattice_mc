add_library(utils_library 
            io_utils.cc
            vector_utils.cc
            ${HEADER_FRUSA_UTILITY})

target_include_directories(utils_library PUBLIC 
                           ${INCLUDE_FRUSA_UTILITY})

target_link_libraries(utils_library PUBLIC compiler_flags)
