### Create the main model library

add_library(model_library model.cc
            ${HEADER_FRUSA_MODELS}
            ${HEADER_FRUSA_UTILITY})

set(MODEL_TYPE
    "fields_chain" CACHE STRING
    "Defines the type of model library to compile"
    PARENT_SCOPE)

target_compile_definitions(model_library PRIVATE "MODEL_TYPE")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fields)

target_include_directories(model_library PUBLIC 
                           ${INCLUDE_FRUSA_MODELS} 
                           ${PROJECT_SOURCE_DIR}/include/models/fields)

target_link_libraries(model_library PUBLIC compiler_flags)
target_link_libraries(model_library PUBLIC fields_library)

