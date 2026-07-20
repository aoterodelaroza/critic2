## CPack configuration for the critic2 installation packages.
set(CPACK_PACKAGE_NAME "critic2")
set(CPACK_PACKAGE_VENDOR "Alberto Otero de la Roza")
set(CPACK_PACKAGE_CONTACT "${critic2_EMAIL}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${critic2_DESCRIPTION}")
set(CPACK_PACKAGE_HOMEPAGE_URL "${critic2_URL}")
set(CPACK_PACKAGE_VERSION "${critic2_VERSION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "critic2")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
set(CPACK_STRIP_FILES ON)

## default generators per platform
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
## CPACK_CREATE_DESKTOP_LINKS names the subset that also get a desktop icon
## (gated on the "Create desktop icon" checkbox); on a GUI build that is both,
## so together with the software-rendering entry below all three land on the
## desktop and in the Start Menu.
if (WIN32 AND ENABLE_GUI)
  set(CPACK_PACKAGE_EXECUTABLES
    "critic2-gui" "critic2"
    "critic2"     "critic2 (console)")
  set(CPACK_CREATE_DESKTOP_LINKS "critic2-gui" "critic2")
elseif (WIN32)
  set(CPACK_PACKAGE_EXECUTABLES "critic2" "critic2")
  set(CPACK_CREATE_DESKTOP_LINKS "critic2")
endif()

## Software-rendering shortcut.
if (WIN32 AND CRITIC2_HAVE_MESA)
  set(CPACK_NSIS_CREATE_ICONS_EXTRA
    "CreateShortCut '$SMPROGRAMS\\$STARTMENU_FOLDER\\critic2 (software rendering).lnk' '$INSTDIR\\critic2 (software rendering).exe'
  StrCmp '$INSTALL_DESKTOP' '1' 0 +2
    CreateShortCut '$DESKTOP\\critic2 (software rendering).lnk' '$INSTDIR\\critic2 (software rendering).exe'")
  set(CPACK_NSIS_DELETE_ICONS_EXTRA
    "Delete '$SMPROGRAMS\\$MUI_TEMP\\critic2 (software rendering).lnk'
  StrCmp '$INSTALL_DESKTOP' '1' 0 +2
    Delete '$DESKTOP\\critic2 (software rendering).lnk'")
endif()

## Every component (program, data files) is mandatory: drop the
## component-selection page altogether.
set(CPACK_MONOLITHIC_INSTALL ON)

include(CPack)
