include( ExternalProject )

ExternalProject_Add( YAML-CPP
  GIT_REPOSITORY "https://github.com/jbeder/yaml-cpp.git"
  GIT_TAG "yaml-cpp-0.6.2"
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  PREFIX ${CMAKE_BINARY_DIR}/external
  INSTALL_DIR ${CMAKE_BINARY_DIR}/external/yaml
  CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

install( DIRECTORY ${CMAKE_BINARY_DIR}/external/yaml/
  DESTINATION . USE_SOURCE_PERMISSIONS )

set( YAML_INCLUDE_DIR
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/external/yaml/include>
  $<INSTALL_INTERFACE:include> )

set( YAMLCPP_LIBRARIES )

foreach( lib yaml-cpp )
  list( APPEND YAMLCPP_LIBRARIES
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/external/yaml/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}>
    $<INSTALL_INTERFACE:lib/${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}> )
endforeach()
