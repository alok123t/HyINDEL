install( PROGRAMS ${CMAKE_SOURCE_DIR}/scripts/plusminus.sh
    DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    RENAME plusminus
    )

install( PROGRAMS ${CMAKE_SOURCE_DIR}/scripts/preProcess.sh
    DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    RENAME preProcess
    )

install( PROGRAMS ${CMAKE_SOURCE_DIR}/scripts/postProcess.sh
    DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    RENAME postProcess
    )

install( PROGRAMS ${CMAKE_SOURCE_DIR}/scripts/assembleInsertions.sh
    DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    RENAME assembleInsertions
    )
    