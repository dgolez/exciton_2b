cmake_minimum_required(VERSION 2.8)
project("exciton")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(COMMON_LIBS "")


# ~~ Add Eigen ~~
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# ~~ Paths and Subdirs ~~
include_directories(${CMAKE_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ~~ Libraries ~~
link_directories(${CMAKE_LIBRARY_PATH})

# ~~ Add HDF5 ~~
if (hdf5)
  #find_library(HDF5_LIB hdf5)
  find_package(HDF5)
  set(HDF5_LIB ${LIBS} ${HDF5_LIBRARIES})
  include_directories(${HDF5_INCLUDE_DIR})
  add_definitions("-DCNTR_USE_HDF5")
endif (hdf5)

# ~~ OpenMP and MPI flags ~~
option(omp "Build with OpenMP support" OFF)
if (omp)
    message(STATUS "Building with OpenMP")
    find_package(OpenMP REQUIRED)
    add_definitions("-DCNTR_USE_OMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

option(mpi "Build with OpenMPI support" OFF)
if (mpi)
    message(STATUS "Building with OpenMPI")
    find_package(MPI REQUIRED)
    add_definitions("-DCNTR_USE_MPI")
    set(CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_COMPILE_FLAGS
        "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    set(CMAKE_CXX_LINK_FLAGS
        "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    include_directories(${MPI_INCLUDE_PATH})
    list(APPEND COMMON_LIBS ${MPI_LIBRARIES})
endif ()


add_subdirectory(program)
