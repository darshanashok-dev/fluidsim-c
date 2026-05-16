cmake_minimum_required(VERSION 3.14)
project(FluidSim C)

# Fix for raylib 5.0 requiring old CMake version
set(CMAKE_POLICY_VERSION_MINIMUM 3.0)
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(CMAKE_C_STANDARD 99)

include(FetchContent)
FetchContent_Declare(
    raylib
    URL https://github.com/raysan5/raylib/archive/refs/heads/master.tar.gz
)
FetchContent_MakeAvailable(raylib)

find_package(OpenMP REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ffast-math -march=native")

add_executable(FluidSim 
    src/main.c 
    src/fluid.c
)

target_include_directories(FluidSim PRIVATE include)
target_link_libraries(FluidSim PRIVATE raylib OpenMP::OpenMP_C)
