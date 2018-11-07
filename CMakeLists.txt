cmake_minimum_required( VERSION 3.1 FATAL_ERROR )

project( indel LANGUAGES CXX )

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

set( THREADS_PREFER_PTHREAD_FLAG ON )
find_package( Threads REQUIRED )

set( EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin" )
set( LIBRARY_OUTPUT_PATH    "${CMAKE_BINARY_DIR}/lib" )

add_subdirectory( "${CMAKE_SOURCE_DIR}/external/bamtools" )
add_subdirectory( "${CMAKE_SOURCE_DIR}/source" )

enable_testing()
add_subdirectory( "${CMAKE_SOURCE_DIR}/test" )
