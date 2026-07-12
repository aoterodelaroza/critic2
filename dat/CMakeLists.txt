install(DIRECTORY lib wfc force_fields DESTINATION "${CMAKE_INSTALL_DATADIR}/critic2")
if (ENABLE_GUI)
  install(DIRECTORY helpdoc shaders assets DESTINATION "${CMAKE_INSTALL_DATADIR}/critic2")
endif()

