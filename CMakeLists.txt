cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(StatisticsTools CXX)

message("Project dir is ${PROJECT_SOURCE_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set(SOURCE_FILES_FIT fit.cxx)
set(EXECUTABLE_FIT fit.exe)

set(SOURCE_FILES_PULLS pulls.cxx)
set(EXECUTABLE_PULLS pulls.exe)

set(SOURCE_FILES_PLOT_PULLS plot_pulls.cxx)
set(EXECUTABLE_PLOT_PULLS plot_pulls.exe)

set(SOURCE_FILES_SCAN scan.cxx)
set(EXECUTABLE_SCAN scan.exe)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++14 COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  add_definitions(-std=c++14)
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()

# CERN ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core Gpad Hist HistFactory MathCore Matrix Minuit Net Physics RIO Rint RooFit RooFitCore RooStats TMVA Thread Tree TreePlayer)

# BOOST
find_package(Boost 1.60 COMPONENTS REQUIRED program_options log log_setup system filesystem date_time thread chrono)

# Includes
set(PROJECT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${YAMLCPP_INCLUDE_DIRS} ${YAMLCPP_INCLUDE_DIR})
include(${ROOT_USE_FILE})
include_directories(${PROJECT_INCLUDE_DIRS})

# ATLAS style
include_directories(${PROJECT_SOURCE_DIR}/atlasstyle-00-03-05/)
file(GLOB sources ${PROJECT_SOURCE_DIR}/atlasstyle-00-03-05/*.C)
file(GLOB headers ${PROJECT_SOURCE_DIR}/atlasstyle-00-03-05/*.h)

# Generate dictionary
ROOT_GENERATE_DICTIONARY(G__StatisticsTools utils.hxx ExtendedModel.hxx ExtendedMinimizer.hxx ${headers} LINKDEF LinkDef.h)

# Shared library with generated dictionary
set(PROJECT_LIBRARIES ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${YAMLCPP_LIBRARIES} ${YAMLCPP_LIBRARY})
add_library(StatisticsTools SHARED utils.hxx ExtendedModel.cxx ExtendedMinimizer.cxx G__StatisticsTools.cxx ${sources})
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
