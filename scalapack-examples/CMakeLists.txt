
cmake_minimum_required(VERSION 3.1...3.21)
project(qrcp LANGUAGES CXX)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")

find_package(BLAS REQUIRED)

find_package(SCALAPACK QUIET)

if (NOT ${SCALAPACK_FOUND})
    message("Scalapack not found")
    set(SCALAPACK_LIBRARIES "/opt/ohpc/pub/libs/gnu8/openmpi3/scalapack/2.0.2/lib/libscalapack.so")
endif()

find_package(MPI REQUIRED COMPONENTS CXX REQUIRED)
message("-- MPI library found in " "${MPI_LIBRARIES}")
message("-- MPI include path is" "${MPI_INCLUDE_PATH}")

add_compile_options(-Wall)
include_directories(${MPI_INCLUDE_PATH})

add_executable(qrcp qrcp-example.cpp)
target_link_libraries(qrcp 
        ${MPI_LIBRARIES} 
        ${SCALAPACK_LIBRARIES} 
        ${BLAS_LIBRARIES}
        )

#mpic++ qrcp-example.cpp -L/opt/ohpc/pub/libs/gnu8/openblas/0.3.7/lib/ -L/opt/ohpc/pub/libs/gnu8/openmpi3/scalapack/2.0.2/lib/  -lscalapack -lopenblas 

