cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(StatisticsTools CXX)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

message("Project dir is ${PROJECT_SOURCE_DIR}")
if(NOT DEFINED ENV{AtlasProject}) # not needed when inside ATLAS environment
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
  set(CMAKE_INSTALL_BINDIR ${PROJECT_SOURCE_DIR}/bin)
  set(CMAKE_INSTALL_LIBDIR ${PROJECT_SOURCE_DIR}/lib)
endif()

set(SOURCE_FILES_BREAKDOWN run_breakdown.cxx)
set(EXECUTABLE_BREAKDOWN run_breakdown.exe)

set(SOURCE_FILES_FIT run_fit.cxx)
set(EXECUTABLE_FIT run_fit.exe)

set(SOURCE_FILES_READ readWS.cxx)
set(EXECUTABLE_READ readWS.exe)

set(SOURCE_FILES_READFIT readWS_fit.cxx)
set(EXECUTABLE_READFIT readWS_fit.exe)

set(SOURCE_FILES_SRFRAC cutPostFitBkg.cxx)
set(EXECUTABLE_SRFRAC cutPostFitBkg.exe)

set(SOURCE_FILES_COMBMYY combMyyPlots.cxx)
set(EXECUTABLE_COMBMYY combMyyPlots.exe)

set(SOURCE_FILES_PULLS run_pulls.cxx)
set(EXECUTABLE_PULLS run_pulls.exe)

set(SOURCE_FILES_SCAN run_scan.cxx)
set(EXECUTABLE_SCAN run_scan.exe)

set(SOURCE_FILES_PLOT_LOGSB plot_logsb.cxx)
set(EXECUTABLE_PLOT_LOGSB plot_logsb.exe)

set(SOURCE_FILES_PLOT_PULLS plot_pulls.cxx)
set(EXECUTABLE_PLOT_PULLS plot_pulls.exe)

set(SOURCE_FILES_PLOT_SCAN plot_scan.cxx)
set(EXECUTABLE_PLOT_SCAN plot_scan.exe)

set(SOURCE_FILES_LIMIT compute_limit_from_scan.cxx)
set(EXECUTABLE_LIMIT compute_limit_from_scan.exe)

#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG(-std=c++14 COMPILER_SUPPORTS_CXX14)
#if(COMPILER_SUPPORTS_CXX14)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
#  add_definitions(-std=c++14)
#else()
#  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
#endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")

# Additional modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules" ${CMAKE_MODULE_PATH})

# CERN ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core Gpad Hist HistFactory MathCore Matrix Minuit Net Physics RIO Rint RooFit RooFitCore RooStats TMVA Thread Tree TreePlayer Graf)

# BOOST
find_package(Boost 1.60 COMPONENTS REQUIRED program_options log log_setup system filesystem date_time thread chrono)

# YAML
include(cmake/Modules/YamlCpp.cmake)

# Eigen
include(cmake/Modules/Eigen.cmake)

# Includes
set(PROJECT_INCLUDE_DIRS include ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${YAML_INCLUDE_DIR})

if(NOT DEFINED ENV{AtlasProject}) # not defined when inside ATLAS environment
  include(${ROOT_USE_FILE})
endif()

include_directories(${PROJECT_INCLUDE_DIRS})

# ATLAS style
include_directories(${PROJECT_SOURCE_DIR}/atlasrootstyle/)
set(atlasrootstyle_sources atlasrootstyle/AtlasLabels.C atlasrootstyle/AtlasStyle.C atlasrootstyle/AtlasUtils.C)
set(atlasrootstyle_headers atlasrootstyle/AtlasLabels.h atlasrootstyle/AtlasStyle.h atlasrootstyle/AtlasUtils.h)

# Generate dictionary
ROOT_GENERATE_DICTIONARY(G__StatisticsTools MODULE utils.hxx RooMultiVarGaussianHighPrecision.h RooTwoSidedCBShape.h ExtendedModel.hxx ExtendedMinimizer.hxx ${atlasrootstyle_headers} LINKDEF include/LinkDef.h)

# Shared library with generated dictionary
set(PROJECT_LIBRARIES ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${YAMLCPP_LIBRARIES})
add_library(StatisticsTools SHARED src/utils.cxx src/RooMultiVarGaussianHighPrecision.cxx src/RooTwoSidedCBShape.cxx src/ExtendedModel.cxx src/ExtendedMinimizer.cxx G__StatisticsTools.cxx ${atlasrootstyle_sources})
target_link_libraries(StatisticsTools ${PROJECT_LIBRARIES})
add_dependencies(StatisticsTools YAML-CPP Eigen)

# Create a main program using the library
add_executable(${EXECUTABLE_FIT} ${SOURCE_FILES_FIT})
target_link_libraries(${EXECUTABLE_FIT} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_READ} ${SOURCE_FILES_READ})
target_link_libraries(${EXECUTABLE_READ} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_READFIT} ${SOURCE_FILES_READFIT})
target_link_libraries(${EXECUTABLE_READFIT} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_SRFRAC} ${SOURCE_FILES_SRFRAC})
target_link_libraries(${EXECUTABLE_SRFRAC} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_COMBMYY} ${SOURCE_FILES_COMBMYY})
target_link_libraries(${EXECUTABLE_COMBMYY} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_PULLS} ${SOURCE_FILES_PULLS})
target_link_libraries(${EXECUTABLE_PULLS} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_PLOT_PULLS} ${SOURCE_FILES_PLOT_PULLS})
target_link_libraries(${EXECUTABLE_PLOT_PULLS} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_SCAN} ${SOURCE_FILES_SCAN})
target_link_libraries(${EXECUTABLE_SCAN} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_PLOT_SCAN} ${SOURCE_FILES_PLOT_SCAN})
target_link_libraries(${EXECUTABLE_PLOT_SCAN} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_BREAKDOWN} ${SOURCE_FILES_BREAKDOWN})
target_link_libraries(${EXECUTABLE_BREAKDOWN} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_PLOT_LOGSB} ${SOURCE_FILES_PLOT_LOGSB})
target_link_libraries(${EXECUTABLE_PLOT_LOGSB} StatisticsTools ${PROJECT_LIBRARIES})

add_executable(${EXECUTABLE_LIMIT} ${SOURCE_FILES_LIMIT})
target_link_libraries(${EXECUTABLE_LIMIT} StatisticsTools ${PROJECT_LIBRARIES})
