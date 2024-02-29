# Create the main model library
add_library(model model.cpp utils.cpp)

# Define the option for the model type 
set(MODEL_TYPE "" CACHE STRING "Select the model to build")

# Select the model library based on the selected option
if(MODEL_TYPE STREQUAL "1d_particles")

  #Add the option to the compiler
  set(PARTICLES_1D_OPTION ON)
  target_compile_definitions(model PRIVATE PARTICLES_1D_OPTION)

  #Create the library and link it to the main model library
  add_library(1d_particles_model STATIC 1d_particles.cpp)
  
  target_link_libraries(model PUBLIC 1d_particles_model)

  target_include_directories(1d_particles_model PRIVATE 
  ${PROJECT_SOURCE_DIR}/include)

elseif(MODEL_TYPE STREQUAL "1d_fields")
  
  #Add the option to the compiler
  set(FIELDS_1D_OPTION ON)
  target_compile_definitions(model PRIVATE FIELDS_1D_OPTION)

  #Create the library and link it to the main model library
  add_library(1d_fields_model STATIC 1d_fields.cpp)

  target_link_libraries(model PUBLIC 1d_fields_model)

  target_include_directories(1d_fields_model PRIVATE 
  ${PROJECT_SOURCE_DIR}/include)

else()

  message( FATAL_ERROR "Please select a valid model!")

endif()

# Set the directory with the header files
target_include_directories(model PRIVATE 
  ${PROJECT_SOURCE_DIR}/include)


