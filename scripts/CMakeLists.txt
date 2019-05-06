install( PROGRAMS ${CMAKE_SOURCE_DIR}/scripts/plusminus.sh
    DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    RENAME plusminus
    )

install( PROGRAMS ${CMAKE_SOURCE_DIR}/scripts/pm-preProcess.sh
    DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    RENAME pm-preProcess
    )

install( PROGRAMS ${CMAKE_SOURCE_DIR}/scripts/pm-postProcess.sh
    DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    RENAME pm-postProcess
    )

install( PROGRAMS ${CMAKE_SOURCE_DIR}/scripts/pm-assembleInsertions.sh
    DESTINATION "${CMAKE_INSTALL_BINDIR}" 
    RENAME pm-assembleInsertions
    )
    