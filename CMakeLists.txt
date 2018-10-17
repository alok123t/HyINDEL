cmake_minimum_required( VERSION 3.1 FATAL_ERROR )

project( indel LANGUAGES CXX )

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

set( EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin" )
set( LIBRARY_OUTPUT_PATH    "${CMAKE_BINARY_DIR}/lib" )

set( SOURCES ${CMAKE_SOURCE_DIR}/source/Cluster.cpp )

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_subdirectory( "${CMAKE_SOURCE_DIR}/external/bamtools" )

include_directories( "${CMAKE_SOURCE_DIR}/external/args" )
include_directories( "${CMAKE_SOURCE_DIR}/external/bamtools/src" )

add_executable( indel ${SOURCES} )
target_link_libraries( indel Threads::Threads )
target_link_libraries( indel BamTools )
