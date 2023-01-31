## For GUI compilation with Windows/MSYS
## Use it as: cmake -C ../cmake/toolchains/msys2.cmake ..

set(CMAKE_GENERATOR "MSYS Makefiles" CACHE INTERNAL "" FORCE)
set(BUILD_RECIPE "MSYS2" CACHE INTERNAL "" FORCE)
