## For GUI compilation with MINGW/wine, use:
## Use it as: cmake -C ../cmake/toolchains/wine.cmake ..

set(CMAKE_GENERATOR "Unix Makefiles" CACHE INTERNAL "" FORCE)
set(BUILD_RECIPE "WINE" CACHE INTERNAL "" FORCE)

