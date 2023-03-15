# Initialize autoconf
AC_PREREQ([2.71])
## AC_INIT([critic2],[m4_esyscmd_s(./tools/version.sh)],[aoterodelaroza@gmail.com])
AC_INIT([critic2],[1.0],[aoterodelaroza@gmail.com])
AC_CONFIG_AUX_DIR([config/gnu])
AC_CONFIG_MACRO_DIR([config/m4])
# These are necessary for older versions of autotools
# m4_include([config/m4/libxc.m4])

# Initialize automake
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wno-portability])

# windows cross-compilation (system-specific, overrides default compiler choice)
AC_ARG_ENABLE([win32],
[  --enable-win32        Use specific compiler and flags for windows cross-compilation (32 bit)],
[case "${enableval}" in
  yes)
    win32=true
    FC=i686-w64-mingw32-gfortran
    F77=i686-w64-mingw32-gfortran
    CC=i686-w64-mingw32-gcc
    CFLAGS+=" -static"
    FCFLAGS+=" -static"
    FFFLAGS+=" -static"
    LDFLAGS+=" -static"
    ;;
  no)  win32=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-win32]) ;;
esac],[win32=false])
AC_ARG_ENABLE([win64],
[  --enable-win64        Use specific compiler and flags for windows cross-compilation (64 bit)],
[case "${enableval}" in
  yes)
    win64=true
    FC=x86_64-w64-mingw32-gfortran
    F77=x86_64-w64-mingw32-gfortran
    CC=x86_64-w64-mingw32-gcc
    CFLAGS+=" -static"
    FCFLAGS+=" -static"
    FFFLAGS+=" -static"
    LDFLAGS+=" -static"
    ;;
  no)  win64=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-win64]) ;;
esac],[win64=false])
AM_CONDITIONAL([WIN], [test x$win32 = xtrue -o x$win64 = xtrue])

# Some fortran compiler flags I'll need
AC_LANG([Fortran])
AC_PROG_FC([gfortran ifort])
AC_PROG_F77([gfortran ifort])
AC_FC_FREEFORM
AC_FC_LINE_LENGTH([unlimited])
AC_OPENMP
FFLAGS+=" $OPENMP_FFLAGS"
FCFLAGS+=" $OPENMP_FCFLAGS"
AC_EXEEXT

# C language
AC_LANG_PUSH([C])
AC_PROG_CC
## AC_CONFIG_HEADERS([config.h])
## AC_CHECK_HEADERS([stdio.h string.h stdlib.h GL/gl3w.h])
AC_LANG_POP

# C++ language
AC_LANG_PUSH([C++])
AC_PROG_CXX
AC_LANG_POP

# check for programs
AC_PROG_RANLIB
AC_PROG_AWK(gawk)
AC_PROG_MAKE_SET
AC_PROG_INSTALL

# Debug compilation, system-specific
AC_ARG_ENABLE(debug,
 [AS_HELP_STRING([--enable-debug],[use debug compiler flags for ifort/gfortran...
  This is a shorthand for me because I can't remember all the flags.
  Not portable at all. Do not use it.
  Do not look at it. Go away, shoo shoo! -- aor, default: no])],,[enable_debug=no])
if test "${enable_debug}" = "yes"; then
   case $FC in
   "gfortran"|"i686-w64-mingw32-gfortran"|"x86_64-w64-mingw32-gfortran"|"gfortran-5"|"gfortran-6"|"gfortran-7"|"gfortran-8"|"gfortran-9"|"gfortran-10")
	FCFLAGS=" -g -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -ffree-line-length-0 -cpp"
	#FFLAGS=" -g -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid -fbacktrace -fdump-core -ffree-line-length-0 -cpp"
	CFLAGS=" -g -Og -fbounds-check -Wall -Wunused-parameter -rdynamic"
	CXXFLAGS=" -g -Og -fbounds-check -Wall -Wunused-parameter -rdynamic"
	LDFLAGS+=" -rdynamic"
	;;
   "ifort")
	FCFLAGS=" -g -C -traceback -debug -fpp -check all -check noarg_temp_created -O0"
	#FFLAGS=" -g -C -traceback -debug -fpp -check all -check noarg_temp_created -O0"
	;;
   esac
fi

# check if we need fallow-argument-mismatch
AX_CHECK_COMPILE_FLAG([-fallow-argument-mismatch], [FCFLAGS="$FCFLAGS -fallow-argument-mismatch"])

# libxc
# ACX_LIBXC

# libcint
AC_ARG_WITH([cint-shared],[AS_HELP_STRING([--with-cint-shared],[Location of the cint shared library.])])
if test "${with_cint_shared}" != ""; then
  LDFLAGS+=" -L${with_cint_shared} -lcint -lblas"
  AC_DEFINE(HAVE_CINT, 1, [Defined if you have the CINT library.])
else
  AC_ARG_WITH([cint-static],[AS_HELP_STRING([--with-cint-static],[Location of the cint static library.])])
  if test "${with_cint_static}" != ""; then
    LIBADD+=" ${with_cint_static}/libcint.a -lblas"
    AC_DEFINE(HAVE_CINT, 1, [Defined if you have the CINT library.])
  fi
fi

# more variables
_os=$(uname -a)
_date=$(date)
AC_SUBST(DEBUG, ${enable_debug})
AC_SUBST(TARGET, ${_os})
AC_SUBST(DATE, ${_date})
AC_SUBST(LIBADD, ${LIBADD})

# Messages
AC_MSG_NOTICE([FC        = ${FC}])
AC_MSG_NOTICE([F77       = ${F77}])
AC_MSG_NOTICE([CC        = ${CC}])
AC_MSG_NOTICE([CXX       = ${CXX}])
AC_MSG_NOTICE([FCFLAGS   = ${FCFLAGS}])
AC_MSG_NOTICE([FFLAGS    = ${FFLAGS}])
AC_MSG_NOTICE([CFLAGS    = ${CFLAGS}])
AC_MSG_NOTICE([CXXFLAGS  = ${CXXFLAGS}])
AC_MSG_NOTICE([LDFLAGS   = ${LDFLAGS}])
AC_MSG_NOTICE([LIBS      = ${LIBS}])
AC_MSG_NOTICE([LIBADD    = ${LIBADD}])

# Final output
AC_CONFIG_FILES([Makefile src/Makefile src/oldlibs/Makefile src/cubpack/Makefile \
	   src/qhull/Makefile src/spglib/Makefile src/json-fortran/Makefile src/emd/Makefile \
	   tools/Makefile src/config@proc.F90 dat/Makefile \
	   dat/wfc/Makefile dat/lib/Makefile])
AC_OUTPUT
