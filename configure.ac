# Initialize autoconf
AC_PREREQ(2.68)
AC_INIT(critic2,m4_esyscmd_s([./tools/version.sh]),[aoterodelaroza@gmail.com,victor@fluor.quimica.uniovi.es,angel@fluor.quimica.uniovi.es])
AC_CONFIG_AUX_DIR([config/gnu])
AC_CONFIG_MACRO_DIR([config/m4])
# These are necessary for older versions of autotools
m4_include([config/m4/ax_f90_module_extension.m4])
m4_include([config/m4/libxc.m4])
m4_include([config/m4/ax_f90_module_flag.m4])

# Initialize automake
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wno-portability])

# Check the fortran compiler
AC_PROG_FC([ifort gfortran])
AC_PROG_F77([ifort gfortran])
AC_PROG_CC

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

# some compiler flags I'll need
AC_PROG_RANLIB
AC_LANG([Fortran])
AC_FC_FREEFORM
AC_FC_LINE_LENGTH([unlimited])
AX_F90_MODULE_EXTENSION
AX_F90_MODULE_FLAG
AC_OPENMP
FFLAGS+=" $OPENMP_FFLAGS"
FCFLAGS+=" $OPENMP_FCFLAGS"
AC_EXEEXT

# check for programs
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
   "gfortran"|"i686-w64-mingw32-gfortran"|"x86_64-w64-mingw32-gfortran"|"gfortran-5"|"gfortran-6")
	FCFLAGS+=" -g -fcheck=all -fbounds-check -Wall -Wunused-parameter -ffpe-trap=invalid -fbacktrace -fdump-core -ffree-line-length-0 -cpp"
	#FFLAGS+=" -g -fcheck=all -fbounds-check -Wall -Wunused-parameter -ffpe-trap=invalid -fbacktrace -fdump-core -ffree-line-length-0 -cpp"
	;;
   "ifort")
        FCFLAGS+=" -g -CU -C -traceback -fpe0 -debug -fpp -check all -check noarg_temp_created -O0"
        #FFLAGS+=" -g -CU -C -traceback -fpe0 -debug -fpp -check all -check noarg_temp_created -O0"
	;;
   esac
   case $CC in
   "gcc"|"x86_64-linux-gnu-gcc-6"|"gcc-5"|"gcc-6"|"cc")
	CFLAGS+=" -g -fbounds-check -Wall -Wunused-parameter -rdynamic"
	;;
   esac
fi

# libxc
ACX_LIBXC

# qhull
AC_ARG_WITH([qhull_inc],
 [AS_HELP_STRING([--with-qhull-inc],[location of the qhull includes.])])
AC_ARG_WITH([qhull_lib],
 [AS_HELP_STRING([--with-qhull-lib],[location of the qhull library.])])

if test "${with_qhull_inc}" != "" && test "${with_qhull_lib}" != ""; then
  CFLAGS+=" -I${with_qhull_inc}"
  LDFLAGS+=" -L${with_qhull_lib} -lqhull"   
else
  ## use the internal qhull    
  CFLAGS+=" -Iqhull/"
  LIBADD+=" qhull/libqhull_critic.a"
fi

# more variables
_os=$(uname -a)
_date=$(date)
AC_SUBST(DEBUG, ${enable_debug})
AC_SUBST(TARGET, ${_os})
AC_SUBST(DATE, ${_date})
AC_SUBST(MODEXT, ${ax_cv_f90_modext})
AC_SUBST(MODINC, ${ax_cv_f90_modflag})
AC_SUBST(LIBADD, ${LIBADD})

# Messages
AC_MSG_NOTICE([FC        = ${FC}])
AC_MSG_NOTICE([F77       = ${F77}])
AC_MSG_NOTICE([CC        = ${CC}])
AC_MSG_NOTICE([FCFLAGS   = ${FCFLAGS}])
AC_MSG_NOTICE([FFLAGS    = ${FFLAGS}])
AC_MSG_NOTICE([CFLAGS    = ${CFLAGS}])
AC_MSG_NOTICE([LDFLAGS   = ${LDFLAGS}])
AC_MSG_NOTICE([LIBS      = ${LIBS}])
AC_MSG_NOTICE([MODINC    = ${MODINC}])
AC_MSG_NOTICE([MODEXT    = ${MODEXT}])
AC_MSG_NOTICE([LIBADD    = ${LIBADD}])

# Final output
AC_OUTPUT([Makefile src/Makefile src/oldlibs/Makefile src/ciftbx/Makefile src/cubpack/Makefile \
           src/qhull/Makefile src/spglib/Makefile tools/Makefile \
           src/config.F90 dat/Makefile doc/Makefile tests/Makefile])

