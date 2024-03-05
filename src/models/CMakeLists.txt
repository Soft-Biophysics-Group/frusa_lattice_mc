### Create the main model library
add_library(model_library model.cc)

target_include_directories(model_library PUBLIC 
                           "${PROJECT_SOURCE_DIR}/include/models"
                           "${PROJECT_SOURCE_DIR}/include/utility/"
                           "${PROJECT_SOURCE_DIR}/include/models/fields"
                           "${PROJECT_SOURCE_DIR}/include/thirdparty")

target_link_libraries(model_library PUBLIC compiler_flags)
target_link_libraries(model_library PUBLIC fields_library)

