set( SOURCES 
    ${CMAKE_SOURCE_DIR}/source/Cluster.cpp
    ${CMAKE_SOURCE_DIR}/source/ArgsParser.cpp
    ${CMAKE_SOURCE_DIR}/source/ArgsParser.hpp
    )

include_directories( "${CMAKE_SOURCE_DIR}/external/args" )
include_directories( "${CMAKE_SOURCE_DIR}/external/bamtools/src" )

add_executable( indel ${SOURCES} )
target_link_libraries( indel Threads::Threads )
target_link_libraries( indel BamTools )
