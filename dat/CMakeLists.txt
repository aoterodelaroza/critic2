install(DIRECTORY lib wfc DESTINATION "${CMAKE_INSTALL_DATADIR}/critic2" COMPONENT data)
if (ENABLE_GUI)
  install(DIRECTORY helpdoc shaders assets DESTINATION "${CMAKE_INSTALL_DATADIR}/critic2" COMPONENT data)
endif()

