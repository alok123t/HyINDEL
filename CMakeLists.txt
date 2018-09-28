cmake_minimum_required( VERSION 2.8 )

project( indel LANGUAGES CXX )

# g++ >= 7.2 or clang++ >= 5.0
if ( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
  if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS "7.2" )
    message( FATAL_ERROR "Requires g++ at least v7.2" )
  endif()
elseif ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
  if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.0" )
    message( FATAL_ERROR "Requires clang++ at least v5.0" )
  endif() 
else()
  message( FATAL_ERROR "Only the following compilers are supported:
    g++ v7.2 or above
    clang++ v5.0 or above" )
endif()

set( CMAKE_CXX_STANDARD 17 )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17" )

set( EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin" )
set( LIBRARY_OUTPUT_PATH    "${CMAKE_BINARY_DIR}/lib" )

add_subdirectory( "${CMAKE_SOURCE_DIR}/external/bamtools" )

include_directories( "${CMAKE_SOURCE_DIR}/external/args" )
include_directories( "${CMAKE_SOURCE_DIR}/external/bamtools/src" )
include_directories( "${CMAKE_SOURCE_DIR}/external/cpp-taskflow" )

set( SOURCES 
    ${CMAKE_SOURCE_DIR}/source/Cluster.cpp )

add_executable( indel ${SOURCES} )
target_link_libraries( indel BamTools )