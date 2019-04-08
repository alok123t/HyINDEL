set( SOURCES 
    ${CMAKE_SOURCE_DIR}/source/Main.cpp
    ${CMAKE_SOURCE_DIR}/source/Parser.cpp
    ${CMAKE_SOURCE_DIR}/source/Parser.hpp
    ${CMAKE_SOURCE_DIR}/source/Detect.cpp
    ${CMAKE_SOURCE_DIR}/source/Detect.hpp
    ${CMAKE_SOURCE_DIR}/source/Misc.hpp
    ${CMAKE_SOURCE_DIR}/source/Writer.cpp
    ${CMAKE_SOURCE_DIR}/source/Writer.hpp
    )

set ( SOURCES_INS
    ${CMAKE_SOURCE_DIR}/source/Insertions.cpp
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