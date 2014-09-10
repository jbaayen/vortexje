cmake_minimum_required(VERSION 2.8.9)
project(Vortexje)

# The version number.
set(Vortexje_VERSION_MAJOR 0)
set(Vortexje_VERSION_MINOR 2)

# Use C++11.
include("cmake/set_cxx_norm.cmake")
set_cxx_norm(CXX_NORM_CXX11)

# Use Eigen3.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3 3.2.2 REQUIRED)

# Base compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O3 -march=native -Wall") # Optimize for this machine, and turn on compiler warnings.
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")             # Optimize for this architecture.
endif()

# Enable PIC for shared library generation.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Use OpenMP.
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
endif()

# Include directories.
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(.)

# Enable testing.
enable_testing()

# Subdirectories.
add_subdirectory(vortexje)
add_subdirectory(tests)
add_subdirectory(examples)
add_subdirectory(doc)

# Install pkg-config file.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vortexje.pc.in ${CMAKE_CURRENT_BINARY_DIR}/vortexje.pc @ONLY)
  
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/vortexje.pc DESTINATION lib/pkgconfig)

# Build a CPack driven installer package
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE  
    "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "${Vortexje_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${Vortexje_VERSION_MINOR}")
include(CPack)
