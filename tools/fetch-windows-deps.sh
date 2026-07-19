#!/bin/bash
# Download the Windows binary dependencies needed to cross-compile the critic2
# GUI package and print the exact commands to build and package it.
#
# Fetched:
#   * GLFW           (always)   -- the GUI windowing/OpenGL library
#   * Mesa           (optional) -- software-OpenGL fallback for machines with
#                                  no GPU / OpenGL < 3.3
#   * OpenBLAS       (optional) -- fast external LAPACK/BLAS
#   * FreeType       (optional) -- higher-quality GUI font rendering. It has no
#                                  MinGW prebuilt binary, so it is cross-built
#                                  from source into a minimal (dependency-free)
#                                  DLL; needs the MinGW toolchain + cmake and is
#                                  skipped gracefully if either is missing.
#   * readline       (optional) -- interactive command-line editing/history.
#                                  GNU readline does not build for MinGW, so
#                                  this cross-builds wineditline (a BSD-licensed
#                                  readline-API-compatible library) and exposes
#                                  it under a readline-compatible prefix. Also
#                                  needs the toolchain + cmake; skipped if absent.
# OpenMP needs no download: it ships with the MinGW toolchain (libgomp), so it
# is simply enabled in the printed configure line.
#
# The other optional libraries critic2 can use (libxc, HDF5, nlopt, libcint,
# tblite) have no MinGW-compatible prebuilt Windows binaries; using them would
# require cross-building each from source, so they are left OFF here.
#
# On success the script prints a ready-to-paste configure/build/package recipe.
# See also the "Windows builds" section of INSTALL.
#
# Usage:
#   tools/fetch-windows-deps.sh [DEST]
#
# DEST is where the dependencies are unpacked (default: ./windows-deps).
# Environment overrides:
#   GLFW_VERSION      GLFW version to fetch       (default 3.4)
#   MESA_VERSION      Mesa version to fetch       (default 24.1.5; 24.x is Win7-safe)
#   OPENBLAS_VERSION  OpenBLAS version to fetch   (default 0.3.34)
#   FREETYPE_VERSION  FreeType version to build   (default 2.13.3)
#   WINEDITLINE_VERSION  wineditline (readline) version to build (default 2.206)
#   FETCH_MESA        set to 0 to skip Mesa       (default 1)
#   FETCH_LAPACK      set to 0 to skip OpenBLAS   (default 1)
#   FETCH_FREETYPE    set to 0 to skip FreeType   (default 1)
#   FETCH_READLINE    set to 0 to skip readline   (default 1)

set -euo pipefail

# repo root (this script lives in <root>/tools), resolved before we cd to DEST
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

GLFW_VERSION="${GLFW_VERSION:-3.4}"
MESA_VERSION="${MESA_VERSION:-24.1.5}"
OPENBLAS_VERSION="${OPENBLAS_VERSION:-0.3.34}"
FREETYPE_VERSION="${FREETYPE_VERSION:-2.13.3}"
WINEDITLINE_VERSION="${WINEDITLINE_VERSION:-2.206}"
FETCH_MESA="${FETCH_MESA:-1}"
FETCH_LAPACK="${FETCH_LAPACK:-1}"
FETCH_FREETYPE="${FETCH_FREETYPE:-1}"
FETCH_READLINE="${FETCH_READLINE:-1}"
DEST="${1:-$(pwd)/windows-deps}"

# required tools
need=(curl unzip)
[ "$FETCH_MESA" = 1 ] && need+=(7z)
for t in "${need[@]}"; do
   if ! command -v "$t" >/dev/null 2>&1; then
      echo "error: '$t' not found (install: curl, unzip, p7zip-full)" >&2
      exit 1
   fi
done

mkdir -p "$DEST"
cd "$DEST"

glfw_dir="glfw-${GLFW_VERSION}.bin.WIN64"
mesa_dir="mesa-${MESA_VERSION}"
openblas_dir="$DEST/openblas-${OPENBLAS_VERSION}-x64"

# --- GLFW ---
if [ ! -d "$glfw_dir" ]; then
   echo ">> downloading GLFW ${GLFW_VERSION} ..." >&2
   curl -fL --progress-bar -o glfw.zip \
      "https://github.com/glfw/glfw/releases/download/${GLFW_VERSION}/glfw-${GLFW_VERSION}.bin.WIN64.zip"
   unzip -q -o glfw.zip
   rm -f glfw.zip
fi

glfw_inc="$DEST/$glfw_dir/include"
glfw_lib="$DEST/$glfw_dir/lib-mingw-w64/libglfw3dll.a"
if [ ! -f "$glfw_inc/GLFW/glfw3.h" ] || [ ! -f "$glfw_lib" ]; then
   echo "error: GLFW headers or import library missing under $glfw_dir" >&2
   exit 1
fi

# --- Mesa (optional) ---
mesa_min=""
if [ "$FETCH_MESA" = 1 ]; then
   if [ ! -d "$mesa_dir" ]; then
      echo ">> downloading Mesa ${MESA_VERSION} ..." >&2
      curl -fL --progress-bar -o mesa.7z \
         "https://github.com/pal1000/mesa-dist-win/releases/download/${MESA_VERSION}/mesa3d-${MESA_VERSION}-release-mingw.7z"
      7z x -o"$mesa_dir" mesa.7z >/dev/null
      rm -f mesa.7z
   fi
   if [ ! -f "$mesa_dir/x64/opengl32.dll" ]; then
      echo "error: Mesa opengl32.dll missing under $mesa_dir/x64" >&2
      exit 1
   fi
   # critic2 only needs the OpenGL (WGL) DLLs; the full x64 dir also carries
   # OpenCL, Vulkan, VA and OSMesa (hundreds of MB) that would bloat the package.
   # Collect just the GL essentials into a minimal directory to use as MESA_DIR.
   mesa_min="$DEST/$mesa_dir/critic2-gl"
   mkdir -p "$mesa_min"
   for d in opengl32.dll libgallium_wgl.dll libglapi.dll dxil.dll libspirv_to_dxil.dll; do
      [ -f "$mesa_dir/x64/$d" ] && cp -f "$mesa_dir/x64/$d" "$mesa_min/"
   done
fi

# --- OpenBLAS (optional, external LAPACK/BLAS) ---
have_lapack=0
if [ "$FETCH_LAPACK" = 1 ]; then
   if [ ! -f "$openblas_dir/lib/libopenblas.dll.a" ]; then
      echo ">> downloading OpenBLAS ${OPENBLAS_VERSION} ..." >&2
      # the x64 zip is the MinGW build (bin/lib/include at the archive root)
      curl -fL --progress-bar -o openblas.zip \
         "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${OPENBLAS_VERSION}/OpenBLAS-${OPENBLAS_VERSION}-x64.zip"
      rm -rf "$openblas_dir"
      unzip -q -o openblas.zip -d "$openblas_dir"
      rm -f openblas.zip
   fi
   if [ ! -f "$openblas_dir/lib/libopenblas.dll.a" ] || [ ! -f "$openblas_dir/bin/libopenblas.dll" ]; then
      echo "error: OpenBLAS import library or DLL missing under $openblas_dir" >&2
      exit 1
   fi
   have_lapack=1
fi

# --- FreeType (optional, cross-built from source) ---
# imgui_freetype only needs FreeType's core API, so we build a minimal DLL with
# zlib/bzip2/png/harfbuzz/brotli all disabled: the result is a single small
# libfreetype DLL with no extra runtime dependencies.
ft_root="$DEST/freetype-${FREETYPE_VERSION}-mingw"
have_freetype=0

build_freetype() {
   local src="$DEST/freetype-${FREETYPE_VERSION}"
   local tar="freetype-${FREETYPE_VERSION}.tar.xz"
   echo ">> downloading and cross-building FreeType ${FREETYPE_VERSION} ..." >&2
   curl -fL --progress-bar -o "$tar" \
      "https://download.savannah.gnu.org/releases/freetype/${tar}" || return 1
   rm -rf "$src"
   tar xf "$tar" || return 1
   rm -f "$tar"
   cmake -S "$src" -B "$src/build-mingw" \
      -DCMAKE_TOOLCHAIN_FILE="$REPO_ROOT/cmake/toolchains/x86_64-w64-mingw32.cmake" \
      -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_INSTALL_PREFIX="$ft_root" \
      -DFT_DISABLE_ZLIB=ON -DFT_DISABLE_BZIP2=ON -DFT_DISABLE_PNG=ON \
      -DFT_DISABLE_HARFBUZZ=ON -DFT_DISABLE_BROTLI=ON \
      >"$DEST/freetype-build.log" 2>&1 || return 1
   cmake --build "$src/build-mingw" -j >>"$DEST/freetype-build.log" 2>&1 || return 1
   cmake --install "$src/build-mingw" >>"$DEST/freetype-build.log" 2>&1 || return 1
   [ -f "$ft_root/lib/libfreetype.dll.a" ] || return 1
}

if [ "$FETCH_FREETYPE" = 1 ]; then
   mingw_cc="$(command -v x86_64-w64-mingw32-gcc-posix 2>/dev/null \
            || command -v x86_64-w64-mingw32-gcc 2>/dev/null || true)"
   if [ -f "$ft_root/lib/libfreetype.dll.a" ]; then
      have_freetype=1                              # already built
   elif [ -z "$mingw_cc" ]; then
      echo ">> skipping FreeType: MinGW-w64 toolchain (x86_64-w64-mingw32-gcc) not found" >&2
   elif ! command -v cmake >/dev/null 2>&1; then
      echo ">> skipping FreeType: cmake not found" >&2
   elif build_freetype; then
      have_freetype=1
   else
      echo "warning: FreeType cross-build failed (see $DEST/freetype-build.log);" >&2
      echo "         continuing without it -- the GUI will use its built-in font rasterizer" >&2
   fi
fi

# --- readline (optional, via wineditline, cross-built from source) ---
# GNU readline needs termcap/termios and does not build for Windows. wineditline
# is a BSD-licensed, readline-API-compatible replacement (its edit.dll exports
# readline/add_history/read_history/write_history, exactly what critic2 links).
# We build it and expose it under a readline-compatible prefix (lib/libreadline
# import lib + include/readline/readline.h) so critic2's FindREADLINE finds it.
rl_root="$DEST/wineditline-${WINEDITLINE_VERSION}-mingw"
have_readline=0

build_readline() {
   local src="$DEST/wineditline-${WINEDITLINE_VERSION}"
   local zip="wineditline-${WINEDITLINE_VERSION}.zip"
   echo ">> downloading and cross-building wineditline ${WINEDITLINE_VERSION} (readline) ..." >&2
   curl -fL --progress-bar -o "$zip" \
      "https://downloads.sourceforge.net/project/mingweditline/${zip}" || return 1
   rm -rf "$src"
   unzip -q -o "$zip" || return 1
   rm -f "$zip"
   # MinGW-on-Linux is case-sensitive: the header is strsafe.h, not Strsafe.h
   sed -i 's/#include <Strsafe.h>/#include <strsafe.h>/' "$src/src/fn_complete.c" || return 1
   cmake -S "$src" -B "$src/build-mingw" \
      -DCMAKE_TOOLCHAIN_FILE="$REPO_ROOT/cmake/toolchains/x86_64-w64-mingw32.cmake" \
      -DCMAKE_BUILD_TYPE=Release >"$DEST/readline-build.log" 2>&1 || return 1
   cmake --build "$src/build-mingw" --target edit -j >>"$DEST/readline-build.log" 2>&1 || return 1
   # wineditline's install() hardcodes bin64/ (pre-populated with MSVC binaries),
   # so assemble a clean readline-compatible prefix from the fresh build output
   rm -rf "$rl_root"
   mkdir -p "$rl_root/bin" "$rl_root/lib" "$rl_root/include/readline"
   cp -f "$src/build-mingw/src/edit.dll"      "$rl_root/bin/"                        || return 1
   cp -f "$src/build-mingw/src/libedit.dll.a" "$rl_root/lib/libreadline.dll.a"       || return 1
   cp -f "$src/src/editline/readline.h"       "$rl_root/include/readline/readline.h" || return 1
   [ -f "$rl_root/lib/libreadline.dll.a" ]
}

if [ "$FETCH_READLINE" = 1 ]; then
   mingw_cc="$(command -v x86_64-w64-mingw32-gcc-posix 2>/dev/null \
            || command -v x86_64-w64-mingw32-gcc 2>/dev/null || true)"
   if [ -f "$rl_root/lib/libreadline.dll.a" ]; then
      have_readline=1                              # already built
   elif [ -z "$mingw_cc" ]; then
      echo ">> skipping readline: MinGW-w64 toolchain (x86_64-w64-mingw32-gcc) not found" >&2
   elif ! command -v cmake >/dev/null 2>&1; then
      echo ">> skipping readline: cmake not found" >&2
   elif build_readline; then
      have_readline=1
   else
      echo "warning: readline (wineditline) cross-build failed (see $DEST/readline-build.log);" >&2
      echo "         continuing without it -- the CLI will use the plain console line editor" >&2
   fi
fi

# --- report: dependencies and the exact build/package recipe ---
# collect the extra find-root prefixes (OpenBLAS, FreeType, readline) into one ';' list
find_roots=()
[ "$have_lapack" = 1 ]   && find_roots+=("$openblas_dir")
[ "$have_freetype" = 1 ] && find_roots+=("$ft_root")
[ "$have_readline" = 1 ] && find_roots+=("$rl_root")
_oldifs=$IFS; IFS=';'; find_roots_joined="${find_roots[*]}"; IFS=$_oldifs

# the still-unavailable optional libraries (built ones drop off this list)
notfetched="libxc, HDF5, nlopt, libcint, tblite"
[ "$have_readline" = 1 ] || notfetched="readline, $notfetched"
[ "$have_freetype" = 1 ] || notfetched="freetype, $notfetched"

{
   echo ""
   echo "======================================================================"
   echo "Windows dependencies ready in: $DEST"
   echo "  GLFW     ${GLFW_VERSION}"
   [ -n "$mesa_min" ]      && echo "  Mesa     ${MESA_VERSION} (software-OpenGL fallback)"
   [ "$have_lapack" = 1 ]  && echo "  OpenBLAS ${OPENBLAS_VERSION} (external LAPACK/BLAS)"
   [ "$have_freetype" = 1 ] && echo "  FreeType ${FREETYPE_VERSION} (GUI font rendering)"
   [ "$have_readline" = 1 ] && echo "  readline ${WINEDITLINE_VERSION} (wineditline; CLI editing/history)"
   echo ""
   echo "Optional libraries NOT included (no MinGW prebuilt binaries -- would"
   echo "need cross-building from source): ${notfetched}."
   echo "They stay OFF in the configure line below."
   echo "======================================================================"
   echo ""
   echo "Run these commands from the critic2 source root to build and package"
   echo "the portable ZIP and the NSIS installer (copy-paste the whole block):"
   echo ""
   echo "  # 1. configure (Release, GUI, software fallback + fast BLAS + OpenMP)"
   echo "  cmake -S . -B build-win \\"
   echo "    -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/x86_64-w64-mingw32.cmake \\"
   echo "    -DCMAKE_BUILD_TYPE=Release \\"
   echo "    -DENABLE_GUI=ON -DENABLE_OPENMP=ON \\"
   echo "    -DGLFW3_INCLUDE_DIR=$glfw_inc \\"
   echo "    -DGLFW3_LIBRARY=$glfw_lib \\"
   if [ -n "$mesa_min" ]; then
      echo "    -DMESA_DIR=$mesa_min \\"
   fi
   if [ "$have_lapack" = 1 ]; then
      echo "    -DUSE_EXTERNAL_LAPACK=ON -DBLA_VENDOR=OpenBLAS \\"
   else
      echo "    -DUSE_EXTERNAL_LAPACK=OFF \\"
   fi
   if [ "$have_freetype" = 1 ]; then
      echo "    -DUSE_FREETYPE=ON \\"
   else
      echo "    -DUSE_FREETYPE=OFF \\"
   fi
   if [ "$have_readline" = 1 ]; then
      echo "    -DUSE_READLINE=ON \\"
   else
      echo "    -DUSE_READLINE=OFF \\"
   fi
   if [ -n "$find_roots_joined" ]; then
      echo "    -DCRITIC2_FIND_ROOT=\"$find_roots_joined\" \\"
   fi
   echo "    -DUSE_LIBXC=OFF -DUSE_HDF5=OFF -DUSE_NLOPT=OFF \\"
   echo "    -DUSE_LIBCINT=OFF -DUSE_TBLITE=OFF"
   echo ""
   echo "  # 2. build (produces bin/critic2.exe and bin/critic2-gui.exe)"
   echo "  cmake --build build-win -j"
   echo ""
   echo "  # 3. package -> build-win/critic2-<version>-win64.{zip,exe}"
   echo "  cmake --build build-win --target package        # both ZIP and NSIS"
   echo "  #   or, from inside build-win/:"
   echo "  #     cpack -G ZIP        # portable"
   echo "  #     cpack -G NSIS       # installer (needs makensis)"
   echo ""
} >&2
