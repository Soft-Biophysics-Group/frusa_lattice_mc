set(HEADER_FRUSA_MODELS 
    ${HEADER_FRUSA_MODELS}
    ${CMAKE_CURRENT_SOURCE_DIR}/fields_state.h
    ${CMAKE_CURRENT_SOURCE_DIR}/fields_interactions.h
    PARENT_SCOPE)

if("${MODEL_TYPE}" STREQUAL "fields_chain")
  set(HEADER_FRUSA_MODELS 
      ${HEADER_FRUSA_MODELS}
      ${CMAKE_CURRENT_SOURCE_DIR}/fields_chain.h
      PARENT_SCOPE)

elseif("${MODEL_TYPE}" STREQUAL "fields_square")

  set(HEADER_FRUSA_MODELS 
      ${HEADER_FRUSA_MODELS}
      ${CMAKE_CURRENT_SOURCE_DIR}/fields_square.h
      PARENT_SCOPE)

elseif("${MODEL_TYPE}" STREQUAL "fields_hexagonal")
  
  set(HEADER_FRUSA_MODELS 
      ${HEADER_FRUSA_MODELS}
      ${CMAKE_CURRENT_SOURCE_DIR}/fields_hexagonal.h
      PARENT_SCOPE)

endif()
