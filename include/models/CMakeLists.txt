#if("${MODEL_TYPE}" STREQUAL "fields_chain" OR
#   "${MODEL_TYPE}" STREQUAL "fields_square" OR
#   "${MODEL_TYPE}" STREQUAL "fields_hexagonal")

# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fields)

#endif()

set(HEADER_FRUSA_MODELS 
    ${HEADER_FRUSA_MODELS}
    ${CMAKE_CURRENT_SOURCE_DIR}/model.h
    PARENT_SCOPE)
