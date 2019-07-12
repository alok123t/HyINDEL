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

include_directories( "${CMAKE_SOURCE_DIR}/external/" )
include_directories( "${CMAKE_SOURCE_DIR}/external/bamtools/src" )

add_executable( HyINDEL_dels ${SOURCES} )
target_link_libraries( HyINDEL_dels Threads::Threads )
target_link_libraries( HyINDEL_dels BamTools )

add_executable( HyINDEL_ins ${SOURCES_INS} )
target_link_libraries( HyINDEL_ins Threads::Threads )

install( TARGETS HyINDEL_dels
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    )

install( TARGETS HyINDEL_ins
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    )