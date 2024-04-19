### CMake Version ##############################################################

cmake_minimum_required(VERSION 3.16...3.22)

################################################################################



### Meta information ###########################################################

set(META_PROJECT_NAME  "frusa_mc")
set(META_DESCRIPTION   "Monte-Carlo framework for models of frustrated self-assembly")
set(META_AUTHOR        "Andrey Zelenskiy")
set(META_AUTHOR_EMAIL  "andrey.zelenskiy@universite-paris-saclay.fr")
set(META_VERSION_MAJOR "1")
set(META_VERSION_MINOR "0")
set(META_VERSION_PATCH "0")
set(META_VERSION       "${META_VERSION_MAJOR}.${META_VERSION_MINOR}.${META_VERSION_PATCH}")

################################################################################



### Create frusa_mc project ####################################################

message( STATUS ">> Building ${META_PROJECT_NAME} executable v${META_VERSION}")

project(${META_PROJECT_NAME} 
        VERSION ${META_VERSION}
        DESCRIPTION ${META_DESCRIPTION}
        LANGUAGES CXX)

################################################################################



### Compiler options ###########################################################

add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_11)
add_compile_options(-Weffc++ -Wextra -Wconversion -Wsign-conversion -Werror)
set (CMAKE_CXX_STANDARD 20)

################################################################################



################################################################################

set(MODEL_TYPE 
    "default" CACHE STRING
    "User-defined option for model library")

################################################################################



### Gather headers and sources #################################################

### Header Groups
set(HEADER_FRUSA_ENGINE)
set(HEADER_FRUSA_MODELS)
set(HEADER_FRUSA_UTILITY)
set(HEADER_FRUSA_THIRDPARTY)
set(HEADER_FRUSA_GEOMETRY)

### Header directories
set(INCLUDE_FRUSA_ENGINE ${PROJECT_SOURCE_DIR}/include/engine)
set(INCLUDE_FRUSA_MODELS ${PROJECT_SOURCE_DIR}/include/models)
set(INCLUDE_FRUSA_UTILITY ${PROJECT_SOURCE_DIR}/include/utility)
set(INCLUDE_FRUSA_THIRDPARTY ${PROJECT_SOURCE_DIR}/include/thirdparty)
set(INCLUDE_FRUSA_GEOMETRY ${PROJECT_SOURCE_DIR}/include/geometry)

### Propagate variables to project subdirectories
add_subdirectory(${PROJECT_SOURCE_DIR}/app)
add_subdirectory(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/engine)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/models)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/utility)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/geometry)

###############################################################################
