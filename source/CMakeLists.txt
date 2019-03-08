set( SOURCES 
    ${CMAKE_SOURCE_DIR}/source/Cluster.cpp
    ${CMAKE_SOURCE_DIR}/source/ArgsParser.cpp
    ${CMAKE_SOURCE_DIR}/source/ArgsParser.hpp
    ${CMAKE_SOURCE_DIR}/source/InputReader.cpp
    ${CMAKE_SOURCE_DIR}/source/InputReader.hpp
    ${CMAKE_SOURCE_DIR}/source/Misc.hpp
    ${CMAKE_SOURCE_DIR}/source/OutputWriter.cpp
    ${CMAKE_SOURCE_DIR}/source/OutputWriter.hpp
    )

set ( SOURCES_INS
    ${CMAKE_SOURCE_DIR}/source/AlignInsertions.cpp
    )

include_directories( "${CMAKE_SOURCE_DIR}/external/args" )
include_directories( "${CMAKE_SOURCE_DIR}/external/bamtools/src" )
include_directories( "${CMAKE_SOURCE_DIR}/external/transwarp/src" )

add_executable( indel ${SOURCES} )
target_link_libraries( indel Threads::Threads )
target_link_libraries( indel BamTools )

add_executable( ins ${SOURCES_INS} )

install( TARGETS indel
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    )

install( TARGETS ins
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    )