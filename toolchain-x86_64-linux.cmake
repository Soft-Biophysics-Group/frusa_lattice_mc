set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR x86_64)

# Set cross compiler
set(CMAKE_C_COMPILER /opt/local/bin/x86_64-elf-gcc)
set(CMAKE_CXX_COMPILER /opt/local/bin/x86_64-elf-g++)
set(CMAKE_ASM_COMPILER /opt/local/bin/x86_64-elf-as)
set(CMAKE_LINKER /opt/local/bin/x86_64-elf-ld)

# Set sysroot if needed
set(CMAKE_SYSROOT /opt/local/x86_64-elf)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

# Prevent CMake from testing the compiler with the wrong assembler
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -D_GLIBCXX_USE_CXX20_FEATURES")
