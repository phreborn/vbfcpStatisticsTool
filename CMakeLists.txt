cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(StatisticsTools CXX)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

message("Project dir is ${PROJECT_SOURCE_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_INSTALL_BINDIR ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_LIBDIR ${PROJECT_SOURCE_DIR}/lib)

set(SOURCE_FILES_FIT fit.cxx)
set(EXECUTABLE_FIT fit.exe)

set(SOURCE_FILES_PULLS pulls.cxx)
set(EXECUTABLE_PULLS pulls.exe)

set(SOURCE_FILES_PLOT_PULLS plot_pulls.cxx)
set(EXECUTABLE_PLOT_PULLS plot_pulls.exe)

set(SOURCE_FILES_SCAN scan.cxx)
set(EXECUTABLE_SCAN scan.exe)

set(SOURCE_FILES_PLOT_SCAN plot_scan.cxx)
set(EXECUTABLE_PLOT_SCAN plot_scan.exe)

set(SOURCE_FILES_BREAKDOWN breakdown.cxx)
set(EXECUTABLE_BREAKDOWN breakdown.exe)

set(SOURCE_FILES_PLOT_LOGSB plot_logsb.cxx)
set(EXECUTABLE_PLOT_LOGSB plot_logsb.exe)

set(SOURCE_FILES_LIMIT compute_limit_from_scan.cxx)
set(EXECUTABLE_LIMIT compute_limit_from_scan.exe)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++14 COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  add_definitions(-std=c++14)
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")

# Additional modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules" ${CMAKE_MODULE_PATH})

# CERN ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core Gpad Hist HistFactory MathCore Matrix Minuit Net Physics RIO Rint RooFit RooFitCore RooStats TMVA Thread Tree TreePlayer)

# BOOST
find_package(Boost 1.60 COMPONENTS REQUIRED program_options log log_setup system filesystem date_time thread chrono)

# YAML
find_package(YamlCpp REQUIRED)

# Eigen
find_package(Eigen 3.0 REQUIRED)
include_directories(${Eigen_INCLUDE_DIRS})

# Includes
set(PROJECT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${YAMLCPP_INCLUDE_DIRS} ${YAMLCPP_INCLUDE_DIR})
include(${ROOT_USE_FILE})
include_directories(${PROJECT_INCLUDE_DIRS})

# ATLAS style
include_directories(${PROJECT_SOURCE_DIR}/atlasstyle-00-03-05/)
file(GLOB sources RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/atlasstyle-00-03-05/*.C)
file(GLOB headers RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/atlasstyle-00-03-05/*.h)

# Generate dictionary
ROOT_GENERATE_DICTIONARY(G__StatisticsTools utils.hxx RooMultiVarGaussianHighPrecision.h ExtendedModel.hxx ExtendedMinimizer.hxx ${headers} LINKDEF LinkDef.h)

# Shared library with generated dictionary
set(PROJECT_LIBRARIES ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${YAMLCPP_LIBRARIES} ${YAMLCPP_LIBRARY})
add_library(StatisticsTools SHARED utils.cxx RooMultiVarGaussianHighPrecision.cxx ExtendedModel.cxx ExtendedMinimizer.cxx G__StatisticsTools.cxx ${sources})
target_link_libraries(StatisticsTools ${PROJECT_LIBRARIES})

# Create a main program using the library
add_executable(${EXECUTABLE_FIT} ${SOURCE_FILES_FIT})
target_link_libraries(${EXECUTABLE_FIT} StatisticsTools ${PROJECT_LIBRARIES})

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
