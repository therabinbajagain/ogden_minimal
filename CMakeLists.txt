cmake_minimum_required(VERSION 3.12)
project(ogden_standalone_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(ogden INTERFACE)
target_include_directories(ogden INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_options(ogden INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
)

add_executable(ogden_3d_uniaxial src/ogden_3d_uniaxial.cpp)
target_link_libraries(ogden_3d_uniaxial PRIVATE ogden)

add_executable(ogden_surface_triangle src/ogden_surface_triangle.cpp)
target_link_libraries(ogden_surface_triangle PRIVATE ogden)

add_executable(test_ogden_material tests/test_ogden_material.cpp)
target_link_libraries(test_ogden_material PRIVATE ogden)

add_executable(test_surface_triangle tests/test_surface_triangle.cpp)
target_link_libraries(test_surface_triangle PRIVATE ogden)

enable_testing()
add_test(NAME ogden_material COMMAND test_ogden_material)
add_test(NAME surface_triangle COMMAND test_surface_triangle)
