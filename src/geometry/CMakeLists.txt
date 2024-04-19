### Collect the header and source files to be compiled
set(SOURCE_GEOMETRY
    ${CMAKE_CURRENT_SOURCE_DIR}/geometry.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/triangular.cc
    )

### Create the geometry library and include the header directories
add_library(geometry_library
            ${SOURCE_GEOMETRY}
            ${HEADER_GEOMETRY}
            ${HEADER_THIRDPARTY}
            ${HEADER_UTILITY})

target_include_directories(geometry_library PUBLIC
                           ${INCLUDE_FRUSA_GEOMETRY}
                           ${INCLUDE_FRUSA_THIRDPARTY}
                           ${INCLUDE_FRUSA_UTILITY})

### Link the library to utils to make sure we have access to all of the
### helper functions
target_link_libraries(geometry_library PRIVATE utils_library)

### Make sure that the compiler flags are the same for all libraries by
### linking them to compiler_flags
target_link_libraries(geometry_library PRIVATE compiler_flags)
