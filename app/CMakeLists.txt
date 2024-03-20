### Create the executable

add_executable(${META_PROJECT_NAME} main.cc
               ${HEADER_FRUSA_MODELS})

### Link the libraries to the executable
target_link_libraries(${META_PROJECT_NAME} PUBLIC compiler_flags)
target_link_libraries(${META_PROJECT_NAME} PUBLIC model_library)
target_link_libraries(${META_PROJECT_NAME} PUBLIC mc_library)

### Setup directories
# Copy python directory
file(COPY ${PROJECT_SOURCE_DIR}/python/
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python)

# Set directory for the input files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/input)
