## CPack configuration for the critic2 installation packages. Included at
## the end of the top-level CMakeLists.txt. Typical use, from the build
## directory of a Release build:
##   cpack            # all default generators for this platform
##   cpack -G TGZ     # a specific generator
## Generators: Linux TGZ (+DEB/RPM when the tools are present), macOS
## TGZ/DragNDrop, Windows ZIP (+NSIS when makensis is present; for
## cross-compilation install the nsis package on the Linux host).

set(CPACK_PACKAGE_NAME "critic2")
set(CPACK_PACKAGE_VENDOR "Alberto Otero de la Roza")
set(CPACK_PACKAGE_CONTACT "${critic2_EMAIL}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${critic2_DESCRIPTION}")
set(CPACK_PACKAGE_HOMEPAGE_URL "${critic2_URL}")
set(CPACK_PACKAGE_VERSION "${critic2_VERSION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "critic2")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
set(CPACK_STRIP_FILES ON)

## default generators per platform, overridable with -DCPACK_GENERATOR=...
if (NOT CPACK_GENERATOR)
  if (WIN32)
    set(CPACK_GENERATOR ZIP)
    find_program(MAKENSIS_EXE makensis)
    mark_as_advanced(MAKENSIS_EXE)
    if (MAKENSIS_EXE)
      list(APPEND CPACK_GENERATOR NSIS)
    endif()
  elseif (APPLE)
    set(CPACK_GENERATOR TGZ DragNDrop)
  else()
    set(CPACK_GENERATOR TGZ)
    find_program(DPKG_SHLIBDEPS_EXE dpkg-shlibdeps)
    mark_as_advanced(DPKG_SHLIBDEPS_EXE)
    if (DPKG_SHLIBDEPS_EXE)
      list(APPEND CPACK_GENERATOR DEB)
    endif()
    find_program(RPMBUILD_EXE rpmbuild)
    mark_as_advanced(RPMBUILD_EXE)
    if (RPMBUILD_EXE)
      list(APPEND CPACK_GENERATOR RPM)
    endif()
  endif()
endif()

## DEB options
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)

## RPM options
set(CPACK_RPM_PACKAGE_LICENSE "GPLv3+")
set(CPACK_RPM_FILE_NAME RPM-DEFAULT)

## NSIS options. CRITIC_HOME is not written to the registry; the binary
## finds its data relative to the executable location.
set(CPACK_NSIS_MODIFY_PATH ON)
set(CPACK_NSIS_DISPLAY_NAME "critic2")
set(CPACK_NSIS_URL_INFO_ABOUT "${critic2_URL}")

## Start Menu / desktop shortcuts. Each pair is the executable name (without
## .exe, found under bin/) and the shortcut label; each shortcut shows the icon
## embedded in its executable. On a GUI build there are two binaries -- the
## windowed critic2-gui.exe (the plain "critic2") and the console critic2.exe
## ("critic2 (console)"). A CLI-only build ships the single console shortcut.
if (WIN32 AND ENABLE_GUI)
  set(CPACK_PACKAGE_EXECUTABLES
    "critic2-gui" "critic2"
    "critic2"     "critic2 (console)")
  set(CPACK_CREATE_DESKTOP_LINKS "critic2-gui")
  ## offer to launch the GUI from the installer's finish page
  set(CPACK_NSIS_MUI_FINISHPAGE_RUN "critic2-gui")
elseif (WIN32)
  set(CPACK_PACKAGE_EXECUTABLES "critic2" "critic2")
endif()

include(CPack)
