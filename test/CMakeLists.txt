set( SOURCES 
    ${CMAKE_SOURCE_DIR}/test/ArgsParser.cpp
    ${CMAKE_SOURCE_DIR}/source/ArgsParser.cpp
    ${CMAKE_SOURCE_DIR}/source/ArgsParser.hpp
    )

include_directories( "${CMAKE_SOURCE_DIR}/source" )
include_directories( "${CMAKE_SOURCE_DIR}/external/args" )

add_executable( TestArgsParser ${SOURCES} )

add_test(NAME ArgsParser COMMAND TestArgsParser )
