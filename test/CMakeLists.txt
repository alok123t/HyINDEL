set( SOURCES 
    ${CMAKE_SOURCE_DIR}/test/Parser.cpp
    ${CMAKE_SOURCE_DIR}/source/Parser.cpp
    ${CMAKE_SOURCE_DIR}/source/Parser.hpp
    )

include_directories( "${CMAKE_SOURCE_DIR}/source" )
include_directories( "${CMAKE_SOURCE_DIR}/external/args" )

add_executable( TestParser ${SOURCES} )
add_test( NAME Parser COMMAND TestParser )
