### Create the main model library

add_library(model_library model.cc
            ${HEADER_FRUSA_MODELS}
            ${HEADER_FRUSA_UTILITY})

target_include_directories(model_library PUBLIC 
                           ${INCLUDE_FRUSA_MODELS})

target_link_libraries(model_library PUBLIC compiler_flags)

### Create a compile option to select model library

set(MODEL_TYPE
    "fields_chain" CACHE STRING
    "Defines the type of model library to compile"
    PARENT_SCOPE)

if("${MODEL_TYPE}" STREQUAL "fields_chain" OR 
   "${MODEL_TYPE}" STREQUAL "fields_square" OR
   "${MODEL_TYPE}" STREQUAL "fields_hexagonal") 

  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fields)

  target_include_directories(model_library PUBLIC 
                             ${PROJECT_SOURCE_DIR}/include/models/fields)

  target_link_libraries(model_library PUBLIC fields_library)
endif()
