#!/bin/bash
# Download the Windows binary dependencies needed to cross-compile the critic2
# GUI package and print the exact commands to build and package it.
#
# Fetched:
#   * GLFW           (always)   -- the GUI windowing/OpenGL library
#   * Mesa           (optional) -- software-OpenGL fallback for machines with
#                                  no GPU / OpenGL < 3.3
#   * OpenBLAS       (optional) -- fast external LAPACK/BLAS
#
# None of the following has a MinGW prebuilt binary, so each is cross-built
# from source. They all need the MinGW toolchain + cmake on the host, and each
# is skipped gracefully (the package is simply built without it) when those are
# missing or the build fails:
#   * FreeType       (optional) -- higher-quality GUI font rendering. Built as a
#                                  minimal, dependency-free DLL.
#   * readline       (optional) -- interactive command-line editing/history.
#                                  GNU readline does not build for MinGW, so
#                                  this cross-builds wineditline (a BSD-licensed
#                                  readline-API-compatible library) and exposes
#                                  it under a readline-compatible prefix.
#   * libxc          (optional) -- exchange-correlation energies and potentials.
#                                  critic2 uses the Fortran 2003 interface, and
#                                  its .mod files have to come from the same
#                                  gfortran that builds critic2 -- so a prebuilt
#                                  binary would not have served here anyway.
#   * nlopt          (optional) -- global optimization: variable-cell structure
#                                  comparison and XRPD profile fitting.
#   * libcint        (optional) -- molecular integrals over Gaussian functions.
#   * xtb            (optional) -- GFN-FF energies and forces for the
#                                  interactive dynamics. critic2 drives it
#                                  through its C API only, so no Fortran module
#                                  files have to match; xtb itself is Fortran
#                                  and pulls in two build-time-only libraries
#                                  of its own (see the xtb section below).
# OpenMP needs no download: it ships with the MinGW toolchain (libgomp), so it
# is simply enabled in the printed configure line.
#
# Two optional libraries are still left OFF. HDF5 decides its type conversions
# by running probe programs, which a cross build cannot execute. tblite does
# cross-build, but it needs a five-library chain of its own (toml-f,
# multicharge, s-dftd3, dftd4 and mctc-lib) and two upstream workarounds -- its
# driver executable cannot link against a DLL that re-exports its own static
# dependencies, and its install rules have no RUNTIME DESTINATION, so the DLL
# never leaves the build tree.
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
#   LIBXC_VERSION     libxc version to build      (default 7.0.0)
#   NLOPT_VERSION     nlopt version to build      (default 2.10.0)
#   LIBCINT_VERSION   libcint version to build    (default 6.1.2)
#   XTB_VERSION       xtb version to build        (default 6.7.1)
#   MCTCLIB_VERSION   mctc-lib version (for xtb)  (default 0.4.1)
#   TESTDRIVE_VERSION test-drive version (for xtb) (default 0.5.0)
#   WINEDITLINE_VERSION  wineditline (readline) version to build (default 2.206)
#   FETCH_MESA        set to 0 to skip Mesa       (default 1)
#   FETCH_LAPACK      set to 0 to skip OpenBLAS   (default 1)
#   FETCH_FREETYPE    set to 0 to skip FreeType   (default 1)
#   FETCH_READLINE    set to 0 to skip readline   (default 1)
#   FETCH_LIBXC       set to 0 to skip libxc      (default 1)
#   FETCH_NLOPT       set to 0 to skip nlopt      (default 1)
#   FETCH_LIBCINT     set to 0 to skip libcint    (default 1)
#   FETCH_XTB         set to 0 to skip xtb        (default 1)
#   JOBS              parallel compile jobs       (default: nproc)

set -euo pipefail

# repo root (this script lives in <root>/tools), resolved before we cd to DEST
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

GLFW_VERSION="${GLFW_VERSION:-3.4}"
MESA_VERSION="${MESA_VERSION:-24.1.5}"
OPENBLAS_VERSION="${OPENBLAS_VERSION:-0.3.34}"
FREETYPE_VERSION="${FREETYPE_VERSION:-2.13.3}"
WINEDITLINE_VERSION="${WINEDITLINE_VERSION:-2.206}"
LIBXC_VERSION="${LIBXC_VERSION:-7.0.0}"
NLOPT_VERSION="${NLOPT_VERSION:-2.10.0}"
LIBCINT_VERSION="${LIBCINT_VERSION:-6.1.2}"
XTB_VERSION="${XTB_VERSION:-6.7.1}"
MCTCLIB_VERSION="${MCTCLIB_VERSION:-0.4.1}"
TESTDRIVE_VERSION="${TESTDRIVE_VERSION:-0.5.0}"
FETCH_MESA="${FETCH_MESA:-1}"
FETCH_LAPACK="${FETCH_LAPACK:-1}"
FETCH_FREETYPE="${FETCH_FREETYPE:-1}"
FETCH_READLINE="${FETCH_READLINE:-1}"
FETCH_LIBXC="${FETCH_LIBXC:-1}"
FETCH_NLOPT="${FETCH_NLOPT:-1}"
FETCH_LIBCINT="${FETCH_LIBCINT:-1}"
FETCH_XTB="${FETCH_XTB:-1}"
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

# --- the cross-built dependencies -----------------------------------------
# Everything from here on has no MinGW prebuilt binary and is built from
# source with the same toolchain that will build critic2. Each is optional:
# if the host lacks the toolchain or cmake, or the build fails, the recipe
# printed at the end simply turns that library OFF.
mingw_cc="$(command -v x86_64-w64-mingw32-gcc-posix 2>/dev/null \
         || command -v x86_64-w64-mingw32-gcc 2>/dev/null || true)"
# Debian splits the MinGW C and Fortran compilers into separate packages, and
# libxc and nlopt both need the Fortran one (ENABLE_FORTRAN / NLOPT_FORTRAN).
# Probed separately so a host with only gcc-mingw-w64 says so up front instead
# of downloading ~70 MB of libxc and failing inside cmake.
mingw_fc="$(command -v x86_64-w64-mingw32-gfortran-posix 2>/dev/null \
         || command -v x86_64-w64-mingw32-gfortran 2>/dev/null || true)"
have_cmake="$(command -v cmake 2>/dev/null || true)"
# Bounded parallelism. "cmake --build -j" with no number passes a bare "-j" to
# make, which is NOT "one job per core" but unlimited: every object of a target
# compiles at once. That is harmless for a small dependency and ruinous for
# libxc, which is one target of ~1000 generated C files.
jobs="${JOBS:-$( (nproc 2>/dev/null || echo 4) )}"

# Configure arguments shared by every cross-built dependency. The policy
# minimum sets a floor, so it is a no-op for any project that already declares
# cmake_minimum_required(VERSION >= 3.5) and applies only to the ones below it
# -- do not narrow it to the dependency that happens to need it today.
# CMake >= 4 dropped compatibility with those older declarations, and without
# this wineditline fails to configure and the package silently loses readline.
cmake_cross=(
   -DCMAKE_TOOLCHAIN_FILE="$REPO_ROOT/cmake/toolchains/x86_64-w64-mingw32.cmake"
   -DCMAKE_BUILD_TYPE=Release
   -DCMAKE_POLICY_VERSION_MINIMUM=3.5
)

# needs_fortran TAG -- true for the dependencies built with a Fortran compiler
needs_fortran() {
   case "$1" in
      libxc|nlopt|xtb) return 0 ;;
      *)           return 1 ;;
   esac
}

# have_markers FILE... -- true when every one of them exists
have_markers() {
   local m
   for m in "$@"; do
      if [ ! -f "$m" ]; then
         return 1
      fi
   done
}

# cross_build TAG NAME LOST MARKER...
#   Builds the dependency by calling build_TAG, which logs to $log. The build
#   is skipped when every MARKER is already there, and counts as failed unless
#   they are all there afterwards -- the same test either way, so a prefix that
#   got as far as the import library but not the Fortran module is never
#   mistaken on the next run for one that is ready to use. Returns 0 when the
#   dependency ends up available; LOST says, in the warning, what the package
#   does without it.
cross_build() {
   local tag="$1" name="$2" lost="$3"
   shift 3
   log="$DEST/$tag-build.log"                     # written by build_$tag
   if have_markers "$@"; then
      return 0                                    # already built
   elif [ -z "$mingw_cc" ]; then
      echo ">> skipping $name: MinGW-w64 toolchain (x86_64-w64-mingw32-gcc) not found" >&2
   elif [ -z "$mingw_fc" ] && needs_fortran "$tag"; then
      echo ">> skipping $name: MinGW-w64 Fortran compiler (x86_64-w64-mingw32-gfortran) not found" >&2
   elif [ -z "$have_cmake" ]; then
      echo ">> skipping $name: cmake not found" >&2
   else
      : >"$log"
      if "build_$tag" && have_markers "$@"; then
         return 0
      fi
      echo "warning: $name cross-build failed (see $log);" >&2
      echo "         continuing without it -- $lost" >&2
   fi
   return 1
}

# fetch_source URL DIR
#   Download and unpack into $DEST, unless DIR is already there. Keeping an
#   existing tree lets an interrupted or failed build resume incrementally
#   (its build-mingw/ lives inside DIR) instead of re-downloading and
#   recompiling from scratch; DIR carries the version, so bumping a version
#   still fetches. Delete DIR to force a clean rebuild.
fetch_source() {
   local url="$1" dir="$2"
   local file="${url##*/}"          # separate: one local expands all its words first
   if [ -d "$dir" ]; then
      return 0
   fi
   # noted in the log so that "see <log>" still explains a download failure,
   # which happens before the build has written anything there
   echo "fetching $url" >>"$log"
   if ! curl -fL --progress-bar -o "$file" "$url"; then
      echo "download failed: $url" >>"$log"
      rm -f "$file"
      return 1
   fi
   case "$file" in
      *.zip) unzip -q -o "$file" || { rm -f "$file"; return 1; } ;;
      *)     tar xf "$file"      || { rm -f "$file"; return 1; } ;;
   esac
   rm -f "$file"
}

# cmake_cross_build SRC PREFIX [-D...]
#   The configure/build/install cycle every dependency but wineditline shares.
cmake_cross_build() {
   local src="$1" prefix="$2"
   shift 2
   cmake -S "$src" -B "$src/build-mingw" "${cmake_cross[@]}" \
      -DCMAKE_INSTALL_PREFIX="$prefix" "$@" >>"$log" 2>&1 || return 1
   cmake --build "$src/build-mingw" -j "$jobs" >>"$log" 2>&1 || return 1
   cmake --install "$src/build-mingw" >>"$log" 2>&1
}

# --- FreeType (optional) ---
# imgui_freetype only needs FreeType's core API, so we build a minimal DLL with
# zlib/bzip2/png/harfbuzz/brotli all disabled: the result is a single small
# libfreetype DLL with no extra runtime dependencies.
ft_root="$DEST/freetype-${FREETYPE_VERSION}-mingw"

build_freetype() {
   local src="$DEST/freetype-${FREETYPE_VERSION}"
   echo ">> downloading and cross-building FreeType ${FREETYPE_VERSION} ..." >&2
   fetch_source "https://download.savannah.gnu.org/releases/freetype/freetype-${FREETYPE_VERSION}.tar.xz" \
      "$src" || return 1
   cmake_cross_build "$src" "$ft_root" \
      -DBUILD_SHARED_LIBS=ON \
      -DFT_DISABLE_ZLIB=ON -DFT_DISABLE_BZIP2=ON -DFT_DISABLE_PNG=ON \
      -DFT_DISABLE_HARFBUZZ=ON -DFT_DISABLE_BROTLI=ON
}

have_freetype=0
if [ "$FETCH_FREETYPE" = 1 ] && cross_build freetype "FreeType" \
   "the GUI will use its built-in font rasterizer" \
   "$ft_root/lib/libfreetype.dll.a"; then
   have_freetype=1
fi

# --- readline (optional, via wineditline) ---
# GNU readline needs termcap/termios and does not build for Windows. wineditline
# is a BSD-licensed, readline-API-compatible replacement (its edit.dll exports
# readline/add_history/read_history/write_history, exactly what critic2 links).
# We build it and expose it under a readline-compatible prefix (lib/libreadline
# import lib + include/readline/readline.h) so critic2's FindREADLINE finds it.
# This one does not use cmake_cross_build: it has no usable install step.
rl_root="$DEST/wineditline-${WINEDITLINE_VERSION}-mingw"

build_readline() {
   local src="$DEST/wineditline-${WINEDITLINE_VERSION}"
   echo ">> downloading and cross-building wineditline ${WINEDITLINE_VERSION} (readline) ..." >&2
   fetch_source "https://downloads.sourceforge.net/project/mingweditline/wineditline-${WINEDITLINE_VERSION}.zip" \
      "$src" || return 1
   # MinGW-on-Linux is case-sensitive: the header is strsafe.h, not Strsafe.h
   sed -i 's/#include <Strsafe.h>/#include <strsafe.h>/' "$src/src/fn_complete.c" || return 1
   cmake -S "$src" -B "$src/build-mingw" "${cmake_cross[@]}" >>"$log" 2>&1 || return 1
   cmake --build "$src/build-mingw" --target edit -j "$jobs" >>"$log" 2>&1 || return 1
   # wineditline's install() hardcodes bin64/ (pre-populated with MSVC binaries),
   # so assemble a clean readline-compatible prefix from the fresh build output
   rm -rf "$rl_root"
   mkdir -p "$rl_root/bin" "$rl_root/lib" "$rl_root/include/readline"
   cp -f "$src/build-mingw/src/edit.dll"      "$rl_root/bin/"                        || return 1
   cp -f "$src/build-mingw/src/libedit.dll.a" "$rl_root/lib/libreadline.dll.a"       || return 1
   cp -f "$src/src/editline/readline.h"       "$rl_root/include/readline/readline.h" || return 1
}

have_readline=0
if [ "$FETCH_READLINE" = 1 ] && cross_build readline "readline (wineditline)" \
   "the CLI will use the plain console line editor" \
   "$rl_root/lib/libreadline.dll.a"; then
   have_readline=1
fi

# --- libxc (optional) ---
# critic2 calls libxc through its Fortran 2003 interface, so the build must
# enable it (ENABLE_FORTRAN) and the resulting xc_f03_lib_m.mod must come from
# the same gfortran that compiles critic2 -- .mod files are not portable
# between compilers or major versions.
libxc_root="$DEST/libxc-${LIBXC_VERSION}-mingw"

build_libxc() {
   local src="$DEST/libxc-${LIBXC_VERSION}"
   echo ">> downloading and cross-building libxc ${LIBXC_VERSION} ..." >&2
   fetch_source "https://gitlab.com/libxc/libxc/-/archive/${LIBXC_VERSION}/libxc-${LIBXC_VERSION}.tar.gz" \
      "$src" || return 1
   cmake_cross_build "$src" "$libxc_root" \
      -DBUILD_SHARED_LIBS=ON -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=ON \
      -DENABLE_FORTRAN=ON -DENABLE_PYTHON=OFF -DBUILD_TESTING=OFF
}

have_libxc=0
if [ "$FETCH_LIBXC" = 1 ] && cross_build libxc "libxc" \
   "the LIBXC keyword will be unavailable" \
   "$libxc_root/lib/libxcf03.dll.a" "$libxc_root/include/xc_f03_lib_m.mod"; then
   have_libxc=1
fi

# --- nlopt (optional) ---
# critic2 uses the F77 interface (include 'nlopt.f'), which is only generated
# and installed when NLOPT_FORTRAN is on -- hence nlopt.f is a marker below,
# not just the import library. The C++ half of nlopt is off: the algorithms
# critic2 asks for (SLSQP, MLSL_LDS, PRAXIS, CCSAQ, MMA) are all in the C core,
# and leaving it out keeps libstdc++ off the CLI's DLL list.
nlopt_root="$DEST/nlopt-${NLOPT_VERSION}-mingw"

build_nlopt() {
   local src="$DEST/nlopt-${NLOPT_VERSION}"
   echo ">> downloading and cross-building nlopt ${NLOPT_VERSION} ..." >&2
   fetch_source "https://github.com/stevengj/nlopt/archive/refs/tags/v${NLOPT_VERSION}.tar.gz" \
      "$src" || return 1
   cmake_cross_build "$src" "$nlopt_root" \
      -DBUILD_SHARED_LIBS=ON \
      -DNLOPT_FORTRAN=ON -DNLOPT_CXX=OFF \
      -DNLOPT_PYTHON=OFF -DNLOPT_OCTAVE=OFF -DNLOPT_MATLAB=OFF \
      -DNLOPT_GUILE=OFF -DNLOPT_SWIG=OFF -DNLOPT_TESTS=OFF
}

have_nlopt=0
if [ "$FETCH_NLOPT" = 1 ] && cross_build nlopt "nlopt" \
   "COMPAREVC and the XRPD profile fits will be unavailable" \
   "$nlopt_root/lib/libnlopt.dll.a" "$nlopt_root/include/nlopt.f"; then
   have_nlopt=1
fi

# --- libcint (optional) ---
# libcint predates C23 and declares old-style function pointers (int (*f)()),
# which gcc >= 15 rejects under its default -std=gnu23, so pin the C dialect.
libcint_root="$DEST/libcint-${LIBCINT_VERSION}-mingw"

build_libcint() {
   local src="$DEST/libcint-${LIBCINT_VERSION}"
   echo ">> downloading and cross-building libcint ${LIBCINT_VERSION} ..." >&2
   fetch_source "https://github.com/sunqm/libcint/archive/refs/tags/v${LIBCINT_VERSION}.tar.gz" \
      "$src" || return 1
   cmake_cross_build "$src" "$libcint_root" \
      -DBUILD_SHARED_LIBS=1 -DCMAKE_C_FLAGS="-std=gnu17" \
      -DWITH_FORTRAN=on -DWITH_CINT2_INTERFACE=on \
      -DENABLE_EXAMPLE=0 -DENABLE_TEST=0
}

have_libcint=0
if [ "$FETCH_LIBCINT" = 1 ] && cross_build libcint "libcint" \
   "molecular integrals over Gaussians will be unavailable" \
   "$libcint_root/lib/libcint.dll.a"; then
   have_libcint=1
fi

# --- xtb (optional) ---
# critic2 drives xtb through its C API only (bind(c,name="xtb_...")), so unlike
# libxc there are no Fortran module files that have to match. xtb itself is
# Fortran, though, and needs two libraries of its own at build time. Those are
# built static and end up absorbed into libxtb.dll, so nothing extra ships.
#
# They have to be supplied rather than left to xtb's own FetchContent fallback:
# that clones mctc-lib's HEAD, which now defines the mctc-lib::mctc-lib alias
# that xtb then tries to create itself, and configuring dies with "another
# target with the same name already exists". Passing *_FIND_METHOD=cmake stops
# the search at the prefix built here, which also pins the versions instead of
# tracking a moving HEAD.
#
# -std=gnu17 is needed for the same reason as libcint: symmetry_i.c defines its
# own true/false, which gcc >= 15 rejects under its default -std=gnu23.
xtb_root="$DEST/xtb-${XTB_VERSION}-mingw"
xtb_deps="$DEST/xtb-deps-mingw"

# xtb_prereq NAME URL DIR -- a static, build-time-only dependency of xtb
xtb_prereq() {
   local name="$1" url="$2" dir="$3"
   if [ -f "$xtb_deps/lib/lib$name.a" ]; then
      return 0
   fi
   echo ">>   ... $name (build-time dependency of xtb)" >&2
   fetch_source "$url" "$dir" || return 1
   cmake_cross_build "$dir" "$xtb_deps" -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF
}

build_xtb() {
   local src="$DEST/xtb-${XTB_VERSION}"
   echo ">> downloading and cross-building xtb ${XTB_VERSION} ..." >&2
   xtb_prereq mctc-lib \
      "https://github.com/grimme-lab/mctc-lib/archive/refs/tags/v${MCTCLIB_VERSION}.tar.gz" \
      "$DEST/mctc-lib-${MCTCLIB_VERSION}" || return 1
   xtb_prereq test-drive \
      "https://github.com/fortran-lang/test-drive/archive/refs/tags/v${TESTDRIVE_VERSION}.tar.gz" \
      "$DEST/test-drive-${TESTDRIVE_VERSION}" || return 1
   fetch_source "https://github.com/grimme-lab/xtb/archive/refs/tags/v${XTB_VERSION}.tar.gz" \
      "$src" || return 1
   cmake_cross_build "$src" "$xtb_root" \
      -DBUILD_SHARED_LIBS=ON -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=ON \
      -DCMAKE_C_FLAGS="-std=gnu17" \
      -DWITH_TBLITE=OFF -DWITH_CPCMX=OFF -DBUILD_TESTING=OFF \
      -DMCTCLIB_FIND_METHOD=cmake -DTESTDRIVE_FIND_METHOD=cmake \
      -DBLA_VENDOR=OpenBLAS \
      -DCRITIC2_FIND_ROOT="$openblas_dir;$xtb_deps"
}

have_xtb=0
if [ "$FETCH_XTB" = 1 ]; then
   if [ "$have_lapack" != 1 ]; then
      # xtb's build does find_package(LAPACK REQUIRED), and in a cross build the
      # only LAPACK around is the OpenBLAS fetched above
      echo ">> skipping xtb: it needs LAPACK/BLAS, and OpenBLAS was not fetched" >&2
   elif cross_build xtb "xtb" "the GFN-FF backend of the dynamics will be unavailable" \
      "$xtb_root/lib/libxtb.dll.a"; then
      have_xtb=1
   fi
fi

# --- report: dependencies and the exact build/package recipe ---
# collect the extra find-root prefixes (OpenBLAS, FreeType, readline) into one ';' list
find_roots=()
[ "$have_lapack" = 1 ]   && find_roots+=("$openblas_dir")
[ "$have_freetype" = 1 ] && find_roots+=("$ft_root")
[ "$have_readline" = 1 ] && find_roots+=("$rl_root")
[ "$have_libxc" = 1 ]    && find_roots+=("$libxc_root")
[ "$have_nlopt" = 1 ]    && find_roots+=("$nlopt_root")
[ "$have_libcint" = 1 ]  && find_roots+=("$libcint_root")
[ "$have_xtb" = 1 ]      && find_roots+=("$xtb_root")
_oldifs=$IFS; IFS=';'; find_roots_joined="${find_roots[*]}"; IFS=$_oldifs

# the still-unavailable optional libraries (built ones drop off this list)
notfetched=""
missing() { notfetched="${notfetched:+$notfetched, }$1"; }
[ "$have_freetype" = 1 ] || missing freetype
[ "$have_readline" = 1 ] || missing readline
[ "$have_libxc" = 1 ]    || missing libxc
[ "$have_nlopt" = 1 ]    || missing nlopt
[ "$have_libcint" = 1 ]  || missing libcint
[ "$have_xtb" = 1 ]      || missing xtb
missing "HDF5, tblite"

# useflag NAME HAVE -- one -DUSE_<NAME>=ON/OFF line of the printed recipe
useflag() {
   if [ "$2" = 1 ]; then
      echo "    -DUSE_$1=ON \\"
   else
      echo "    -DUSE_$1=OFF \\"
   fi
}

{
   echo ""
   echo "======================================================================"
   echo "Windows dependencies ready in: $DEST"
   echo "  GLFW     ${GLFW_VERSION}"
   [ -n "$mesa_min" ]      && echo "  Mesa     ${MESA_VERSION} (software-OpenGL fallback)"
   [ "$have_lapack" = 1 ]  && echo "  OpenBLAS ${OPENBLAS_VERSION} (external LAPACK/BLAS)"
   [ "$have_freetype" = 1 ] && echo "  FreeType ${FREETYPE_VERSION} (GUI font rendering)"
   [ "$have_readline" = 1 ] && echo "  readline ${WINEDITLINE_VERSION} (wineditline; CLI editing/history)"
   [ "$have_libxc" = 1 ]    && echo "  libxc    ${LIBXC_VERSION} (exchange-correlation functionals)"
   [ "$have_nlopt" = 1 ]    && echo "  nlopt    ${NLOPT_VERSION} (COMPAREVC, XRPD profile fitting)"
   [ "$have_libcint" = 1 ]  && echo "  libcint  ${LIBCINT_VERSION} (molecular integrals over Gaussians)"
   [ "$have_xtb" = 1 ]      && echo "  xtb      ${XTB_VERSION} (GFN-FF energies and forces)"
   echo ""
   echo "Optional libraries NOT included: ${notfetched}."
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
   useflag FREETYPE "$have_freetype"
   useflag READLINE "$have_readline"
   useflag LIBXC    "$have_libxc"
   useflag NLOPT    "$have_nlopt"
   useflag LIBCINT  "$have_libcint"
   useflag XTB      "$have_xtb"
   if [ -n "$find_roots_joined" ]; then
      echo "    -DCRITIC2_FIND_ROOT=\"$find_roots_joined\" \\"
   fi
   echo "    -DUSE_HDF5=OFF -DUSE_TBLITE=OFF"
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
