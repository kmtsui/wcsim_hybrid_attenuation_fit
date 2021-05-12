cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

#Use the compilers found in the path
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project(OpticalFit)
SET(OPTICALFIT_VERSION "0.0")

#Set this to TRUE to enable build debugging messages
set(BUILD_DEBUG_MSGS TRUE)
include(${CMAKE_SOURCE_DIR}/cmake/cmessage.cmake)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Debug)
elseif(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

cmessage(STATUS "CMAKE_INSTALL_PREFIX: \"${CMAKE_INSTALL_PREFIX}\"")
cmessage(STATUS "CMAKE_BUILD_TYPE: \"${CMAKE_BUILD_TYPE}\"")

# Make a dummy link so that headers can be included sensibly
execute_process(COMMAND ln -s ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/OPTICALFIT OUTPUT_QUIET ERROR_QUIET)

include_directories(${CMAKE_BINARY_DIR})

################################################################################
#                            Check Dependencies
################################################################################

##################################  ROOT  ######################################
include(${CMAKE_SOURCE_DIR}/cmake/ROOTSetup.cmake)

################################## COMPILER ####################################
include(${CMAKE_SOURCE_DIR}/cmake/c++CompilerSetup.cmake)

##################################  WCSIM  ######################################
include(${CMAKE_SOURCE_DIR}/cmake/WCSIMSetup.cmake)

##################################  ####  ######################################

add_subdirectory(src)
add_subdirectory(app)

configure_file(cmake/setup.sh.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" DESTINATION
  ${CMAKE_INSTALL_PREFIX})

include(${CMAKE_SOURCE_DIR}/cmake/opticalfit-config.cmake)