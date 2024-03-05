### Create the executable
add_executable(${META_PROJECT_NAME} main.cpp)

### Link the libraries to the executable
target_link_libraries(${META_PROJECT_NAME} PUBLIC compiler_flags)
target_link_libraries(${META_PROJECT_NAME} PUBLIC model)
#target_link_libraries(${META_PROJECT_NAME} PUBLIC mc_engine)

target_include_directories(${META_PROJECT_NAME} PUBLIC 
  "${PROJECT_SOURCE_DIR}/include/models")

# Add the model library
#add_subdirectory(src/models)

# Set directory with the input files
file(COPY ${PROJECT_SOURCE_DIR}/input/model_params.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/app)

file(COPY ${PROJECT_SOURCE_DIR}/input/mc_params.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

