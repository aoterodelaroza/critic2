## Copyright (C) 2010-2015 M. Marques, X. Andrade, D. Strubbe, M. Oliveira
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id: libxc.m4 13622 2015-03-28 01:31:52Z dstrubbe $
##

AC_DEFUN([ACX_LIBXC], [
acx_libxc_ok=no

dnl Check if the library was given in the command line
dnl if not, use environment variables or defaults
AC_ARG_WITH(libxc-prefix, [AS_HELP_STRING([--with-libxc-prefix=DIR], [Directory where libxc was installed.])])

# Set FCFLAGS_LIBXC only if not set from environment
if test x"$FCFLAGS_LIBXC" = x; then
  case $with_libxc_prefix in
    "") FCFLAGS_LIBXC="$ax_cv_f90_modflag/usr/include" ;;
    *)  FCFLAGS_LIBXC="$ax_cv_f90_modflag$with_libxc_prefix/include" ;;
  esac
fi

AC_ARG_WITH(libxc-include, [AS_HELP_STRING([--with-libxc-include=DIR], [Directory where libxc Fortran headers were installed.])])
case $with_libxc_include in
  "") ;;
  *)  FCFLAGS_LIBXC="$ax_cv_f90_modflag$with_libxc_include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_libxc_save_LIBS="$LIBS"
acx_libxc_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for libxc])

testprog="AC_LANG_PROGRAM([],[
  use xc_f90_lib_m
  implicit none
  integer :: major
  integer :: minor
  call xc_f90_version(major, minor)])"

FCFLAGS="$FCFLAGS_LIBXC $acx_libxc_save_FCFLAGS"

# set from environment variable, if not blank
if test ! -z "$LIBS_LIBXC"; then
  LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
  AC_LINK_IFELSE($testprog, [acx_libxc_ok=yes], [])
fi

if test ! -z "$with_libxc_prefix"; then
  # static linkage, separate Fortran interface (libxc 2.2.0 and later)
  if test x"$acx_libxc_ok" = xno; then
    LIBS_LIBXC="$with_libxc_prefix/lib/libxcf90.a $with_libxc_prefix/lib/libxc.a"
    LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
    AC_LINK_IFELSE($testprog, [acx_libxc_ok=yes], [])
  fi
  
  # static linkage, combined Fortran interface (libxc 2.0.x, 2.1.x)
  if test x"$acx_libxc_ok" = xno; then
    LIBS_LIBXC="$with_libxc_prefix/lib/libxc.a"
    LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
    AC_LINK_IFELSE($testprog, [acx_libxc_ok=yes], [])
  fi
fi

# dynamic linkage, separate Fortran interface (libxc 2.2.0 and later)
if test x"$acx_libxc_ok" = xno; then
  if test ! -z "$with_libxc_prefix"; then
    LIBS_LIBXC="-L$with_libxc_prefix/lib"
  else
    LIBS_LIBXC=""
  fi
  LIBS_LIBXC="$LIBS_LIBXC -lxcf90 -lxc"
  LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
  AC_LINK_IFELSE($testprog, [acx_libxc_ok=yes], [])
fi

# dynamic linkage, combined Fortran interface (libxc 2.0.x, 2.1.x)
if test x"$acx_libxc_ok" = xno; then
  if test ! -z "$with_libxc_prefix"; then
    LIBS_LIBXC="-L$with_libxc_prefix/lib"
  else
    LIBS_LIBXC=""
  fi
  LIBS_LIBXC="$LIBS_LIBXC -lxc"
  LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
  AC_LINK_IFELSE($testprog, [acx_libxc_ok=yes], [])
fi

AC_MSG_RESULT([$acx_libxc_ok ($FCFLAGS_LIBXC $LIBS_LIBXC)])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_libxc_ok" = xyes; then
  AC_DEFINE(HAVE_LIBXC, 1, [Defined if you have the LIBXC library.])
fi

acx_libxc_hyb_mgga_ok=no
AC_MSG_CHECKING([whether libxc has support for hybrid meta-GGAs (>= v 2.1)])

testprog="AC_LANG_PROGRAM([],[
  use xc_f90_lib_m
  implicit none
  integer :: x = XC_FAMILY_HYB_MGGA])"
AC_LINK_IFELSE($testprog, [acx_libxc_hyb_mgga_ok=yes], [])

AC_MSG_RESULT([$acx_libxc_hyb_mgga_ok])
if test x"$acx_libxc_hyb_mgga_ok" = xyes; then
  AC_DEFINE(HAVE_LIBXC_HYB_MGGA, 1, [Defined if LIBXC library has support for hybrid meta-GGAs.])
fi

AC_SUBST(FCFLAGS_LIBXC)
AC_SUBST(LIBS_LIBXC)
if test x"$acx_libxc_ok" = xno; then
  LIBS="$acx_libxc_save_LIBS"
  FCFLAGS="$acx_libxc_save_FCFLAGS"
fi
])dnl ACX_LIBXC
