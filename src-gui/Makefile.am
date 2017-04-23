# binary and sources
bin_PROGRAMS = gcritic2
gcritic2_SOURCES = main.cpp

# subdirectories
SUBDIRS=

# static libraries
gcritic2_LDADD=$(top_srcdir)/src/libcritic2.a \
	$(top_srcdir)/src/oldlibs/libmisc.a \
	$(top_srcdir)/src/ciftbx/libciftbx.a \
	$(top_srcdir)/src/cubpack/libcubpack.a \
	$(top_srcdir)/src/spglib/libspglib.a \
	$(top_srcdir)/src/qhull/libqhull_critic.a

gcritic2_DEPENDENCIES=$(gcritic2_LDADD)

# Include the header files
AM_CXXFLAGS= -I$(top_srcdir)/include
AM_LDFLAGS=-lgfortran -lgomp

