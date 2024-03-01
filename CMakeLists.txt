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

################################################################################



### Gather headers and sources #################################################

### Header groups
set(HEADER_FRUSA_ROOT)
set(HEADER_FRUSA_ENGINE)
set(HEADER_FRUSA_MODELS)
set(HEADER_FRUSA_UTILITY)

### Source groups
set(SOURCE_FRUSA_ROOT)
set(SOURCE_FRUSA_ENGINE)
set(SOURCE_FRUSA_MODELS)
set(SOURCE_FRUSA_UTILITY)

### App group
set(APP_FRUSA)

### Add subdirectories
#add_subdirectory(include)
#add_subdirectory(src)
add_subdirectory(app)

###############################################################################
