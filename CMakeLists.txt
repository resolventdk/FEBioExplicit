cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project( FEBioExplicit )

# Set a default build type if none was specified
set(default_build_type "Release")

add_compile_options(-fPIC)

# FEBIO
if(NOT EXISTS ${FEBIO_LIB_DIR}/libfecore.a)
  message(SEND_ERROR "Could not find FEBio library (libfecore.a). Check FEBIO_LIB_DIR.")
  return()
endif()
ADD_LIBRARY(fecore STATIC IMPORTED)
SET_TARGET_PROPERTIES(fecore PROPERTIES IMPORTED_LOCATION ${FEBIO_LIB_DIR}/libfecore.a)

# if(NOT EXISTS ${FEBIO_LIB_DIR}/libnumcore.a)
#   message(SEND_ERROR "Could not find FEBio library (libnumcore.a). Check FEBIO_LIB_DIR.")
#   return()
# endif()
# ADD_LIBRARY(numcore STATIC IMPORTED)
# SET_TARGET_PROPERTIES(numcore PROPERTIES IMPORTED_LOCATION ${FEBIO_LIB_DIR}/libnumcore.a)

if(NOT EXISTS ${FEBIO_LIB_DIR}/libfebiomech.a)
  message(SEND_ERROR "Could not find FEBio library (libfebiomech.a). Check FEBIO_LIB_DIR.")
  return()
endif()
ADD_LIBRARY(febiomech STATIC IMPORTED)
SET_TARGET_PROPERTIES(febiomech PROPERTIES IMPORTED_LOCATION ${FEBIO_LIB_DIR}/libfebiomech.a)


# OPENMP
find_package(OpenMP)
if ( NOT OPENMP_FOUND )
  message(STATUS "This project requires OpenMP, and will not be compiled.")
  return()  
endif()
add_compile_options(-fopenmp)

# The plugin library
file(GLOB HDR_FEBioExplicit "src/*.h")
file(GLOB SRC_FEBioExplicit "src/*.cpp")
add_library(FEBioExplicit SHARED ${HDR_FEBioExplicit} ${SRC_FEBioExplicit})

target_link_libraries(FEBioExplicit febiomech fecore)
target_include_directories(FEBioExplicit PRIVATE ${FEBIO_ROOT})


