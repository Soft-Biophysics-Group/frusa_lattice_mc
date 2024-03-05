### Create the main model library
add_library(model model.cpp)

target_include_directories(model PUBLIC 
                           "${PROJECT_SOURCE_DIR}/include/models"
                           "${PROJECT_SOURCE_DIR}/include/thirdparty")

target_link_libraries(model PUBLIC compiler_flags)
