cmake_minimum_required(VERSION 3.13)
project(Thesis LANGUAGES CXX CUDA)

# Explicitly enable C++17 and CUDA standards (if needed)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Specify CUDA architectures (adjust based on your GPU)
set(CMAKE_CUDA_ARCHITECTURES 75 80 86)

# Add optimalisation flags
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
# Add executable with both C++ and CUDA sources

add_executable(kernel_tests src/sesion2/main_3.cpp src/sesion2/hyrro_kernel_fixed_2.cu)
add_executable(GPU_VS_CPU src/sesion3/main.cpp src/sesion3/hyrro_kernel.cu src/sesion3/edit_hyrro.cpp)