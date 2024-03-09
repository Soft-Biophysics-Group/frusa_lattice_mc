set(HEADER_FIELDS
    ${INCLUDE_FRUSA_MODELS}/fields/fields_state.h
    ${INCLUDE_FRUSA_MODELS}/fields/fields_interactions.h)

set(SOURCE_FIELDS 
    fields_state.cc
    fields_interactions.cc)

if("${MODEL_TYPE}" STREQUAL "fields_chain")

  set(HEADER_FIELDS 
      ${HEADER_FIELDS} 
      ${INCLUDE_FRUSA_MODELS}/fields/fields_chain.h)
  
  set(SOURCE_FIELDS 
      ${SOURCE_FIELDS} 
      fields_chain.cc)
  
    target_compile_definitions(model_library PUBLIC "FIELDS_CHAIN")

elseif("${MODEL_TYPE}" STREQUAL "fields_square")

  set(HEADER_FIELDS 
      ${HEADER_FIELDS} 
      ${INCLUDE_FRUSA_MODELS}/fields/fields_square.h)
  
  
  set(SOURCE_FIELDS 
      ${SOURCE_FIELDS} 
      fields_square.cc)

    target_compile_definitions(model_library PUBLIC "FIELDS_SQUARE")

elseif("${MODEL_TYPE}" STREQUAL "fields_hexagonal") 

  set(HEADER_FIELDS 
      ${HEADER_FIELDS} 
      ${INCLUDE_FRUSA_MODELS}/fields/fields_hexagonal.h)
  
  set(SOURCE_FIELDS 
      ${SOURCE_FIELDS} 
      fields_hexagonal.cc)
  
    target_compile_definitions(model_library PUBLIC "FIELDS_HEXAGONAL")

endif()

add_library(fields_library ${SOURCE_FIELDS} ${HEADER_FIELDS})

target_include_directories(fields_library PUBLIC
                           ${PROJECT_SOURCE_DIR}/include/models/fields
                           ${INCLUDE_FRUSA_UTILITY}
                           ${INCLUDE_FRUSA_THIRDPARTY})

target_link_libraries(fields_library PUBLIC compiler_flags)
target_link_libraries(fields_library PRIVATE vector_utils)
