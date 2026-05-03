cmake_minimum_required(VERSION 3.10)
project(RDKitExample)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(CMAKE_CXX_STANDARD 20)

# Find RDKit (ensure CONDA_PREFIX is in your path or passed as CMAKE_PREFIX_PATH)
find_package(RDKit REQUIRED)
find_package(Boost REQUIRED)

add_executable(test_interface src/TestInterface.cpp)
target_link_libraries(test_interface -static-libgcc -static-libstdc++ RDKit::GraphMol RDKit::SmilesParse)
target_include_directories(test_interface PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../src")
