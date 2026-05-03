cmake_minimum_required(VERSION 3.10)
project(RDKitExample)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(CMAKE_CXX_STANDARD 20)

# Find RDKit (ensure CONDA_PREFIX is in your path or passed as CMAKE_PREFIX_PATH)
find_package(RDKit REQUIRED)
find_package(Boost REQUIRED)

add_executable(umf-sim-query src/SimilarityQuery.cpp)
target_link_libraries(umf-sim-query -static-libgcc -static-libstdc++ RDKit::GraphMol RDKit::SmilesParse RDKit::Fingerprints RDKit::DataStructs RDKit::GraphMol)
target_include_directories(umf-sim-query PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_options(umf-sim-query PRIVATE -O3)
