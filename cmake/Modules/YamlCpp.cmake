include(ExternalProject)
ExternalProject_Add(YAML-CPP
	        GIT_REPOSITORY  "https://github.com/jbeder/yaml-cpp.git"
	        GIT_TAG  "yaml-cpp-0.6.2"
	        UPDATE_COMMAND ""
	        PATCH_COMMAND ""
	    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/yaml-cpp
	    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}
)

set( YAMLCPP_LIBRARIES )
foreach( lib yaml-cpp )
   list( APPEND YAMLCPP_LIBRARIES
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}>
      $<INSTALL_INTERFACE:lib/${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}> )
endforeach()
