cmake_minimum_required( VERSION 3.0 FATAL_ERROR )

project( HyINDEL LANGUAGES CXX )

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

set( THREADS_PREFER_PTHREAD_FLAG ON )
find_package( Threads REQUIRED )

add_subdirectory( "${CMAKE_SOURCE_DIR}/external/bamtools" )
add_subdirectory( "${CMAKE_SOURCE_DIR}/scripts" )
add_subdirectory( "${CMAKE_SOURCE_DIR}/source" )
