set(HEADER_FIELDS
    ${INCLUDE_FRUSA_MODELS}/fields/fields_parameters.h
    ${INCLUDE_FRUSA_MODELS}/fields/fields_state.h
    ${INCLUDE_FRUSA_MODELS}/fields/fields_geometry.h
    ${INCLUDE_FRUSA_MODELS}/fields/fields_interactions.h
    ${INCLUDE_FRUSA_MODELS}/fields/fields_update.h
    ${INCLUDE_FRUSA_THIRDPARTY}/json.hpp)

set(SOURCE_FIELDS 
    fields_parameters.cc
    fields_state.cc
    fields_geometry.cc
    fields_interactions.cc
    fields_update.cc)

add_library(fields_library 
            ${SOURCE_FIELDS}
            ${HEADER_FIELDS})

if("${MODEL_TYPE}" STREQUAL "fields_chain")

  target_compile_definitions(fields_library PRIVATE "FIELDS_CHAIN")
  
  add_library(user_library STATIC fields_chain.cc)

elseif("${MODEL_TYPE}" STREQUAL "fields_square")

  target_compile_definitions(fields_library PRIVATE "FIELDS_SQUARE")
  
  add_library(user_library STATIC fields_square.cc)

elseif("${MODEL_TYPE}" STREQUAL "fields_hexagonal") 
    
  target_compile_definitions(fields_library PRIVATE "FIELDS_HEXAGONAL")
  
  add_library(user_library STATIC fields_hexagonal.cc)

endif()
  
target_include_directories(user_library PUBLIC
                           ${INCLUDE_FRUSA_MODELS}/fields
                           ${INCLUDE_FRUSA_UTILITY})

target_link_libraries(user_library PRIVATE compiler_flags)
target_link_libraries(user_library PRIVATE vector_utils)

target_link_libraries(fields_library PRIVATE user_library) 

target_include_directories(fields_library PUBLIC
                           ${INCLUDE_FRUSA_MODELS}/fields
                           ${INCLUDE_FRUSA_UTILITY}
                           ${INCLUDE_FRUSA_THIRDPARTY})

target_link_libraries(fields_library PRIVATE compiler_flags)
target_link_libraries(fields_library PRIVATE vector_utils)

