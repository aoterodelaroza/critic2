# binary and sources
bin_PROGRAMS = gcritic2
gcritic2_SOURCES = main.cpp imgui_impl_glfw_gl3.cpp imgui.cpp imgui_draw.cpp \
	matrix_math.cpp gl3w.c geometry.cpp imguifilesystem.cpp \
	guiapps.cpp callback.cpp shader.cpp

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
AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include `pkg-config --cflags glfw3` -Wall -Wformat -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11 -w -g
AM_LDFLAGS=-ldl -lgfortran -lgomp -lglfw -lGL

