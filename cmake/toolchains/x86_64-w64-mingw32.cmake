## Toolchain file for cross-compiling critic2 to 64-bit Windows with
## MinGW-w64 (e.g. Debian/Ubuntu packages gfortran-mingw-w64-x86-64,
## g++-mingw-w64-x86-64). Native MSYS2 builds do NOT need this file;
## configure with -G "MSYS Makefiles" in the MSYS2 shell instead.
##
## Usage:
##   cmake -B build-win -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/x86_64-w64-mingw32.cmake \
##     [-DENABLE_GUI=ON -DGLFW3_ROOT=<path to unpacked glfw-*.bin.WIN64>] \
##     [-DUSE_EXTERNAL_LAPACK=OFF] [-DUSE_LIBXC=OFF] [-DUSE_READLINE=OFF] \
##     [-DCMAKE_CROSSCOMPILING_EMULATOR=wine] \
##     [-DCRITIC2_ROOT=<Windows path to the installed share/ directory>] \
##     [-DCMAKE_INSTALL_PREFIX=<staging directory>]
##
## Libraries are searched only under CMAKE_FIND_ROOT_PATH; pass additional
## prefixes per package (e.g. GLFW3_ROOT) or extend CMAKE_FIND_ROOT_PATH.

set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_SYSTEM_PROCESSOR x86_64)
set(TOOLCHAIN_PREFIX x86_64-w64-mingw32)

set(CMAKE_Fortran_COMPILER ${TOOLCHAIN_PREFIX}-gfortran)
set(CMAKE_C_COMPILER ${TOOLCHAIN_PREFIX}-gcc)
set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)

set(CMAKE_FIND_ROOT_PATH /usr/${TOOLCHAIN_PREFIX})
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
