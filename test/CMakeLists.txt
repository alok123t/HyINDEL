set( SOURCES 
    ${CMAKE_SOURCE_DIR}/test/Parser.cpp
    ${CMAKE_SOURCE_DIR}/source/Parser.cpp
    ${CMAKE_SOURCE_DIR}/source/Parser.hpp
    )

include_directories( "${CMAKE_SOURCE_DIR}/external" )
include_directories( "${CMAKE_SOURCE_DIR}/source" )

add_executable( pm-test ${SOURCES} )
add_test( NAME Parser COMMAND pm-test )
