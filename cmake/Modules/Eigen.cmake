include( ExternalProject )

ExternalProject_Add( Eigen
  GIT_REPOSITORY  "https://github.com/eigenteam/eigen-git-mirror.git"
  GIT_TAG  "3.3.7"
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  PREFIX ${CMAKE_BINARY_DIR}/external
  INSTALL_DIR ${CMAKE_BINARY_DIR}/external/Eigen
  CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

install( DIRECTORY ${CMAKE_BINARY_DIR}/external/Eigen/
  DESTINATION . USE_SOURCE_PERMISSIONS )

set( EIGEN3_INCLUDE_DIR
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/external/Eigen/include/eigen3>
  $<INSTALL_INTERFACE:include/eigen3> )