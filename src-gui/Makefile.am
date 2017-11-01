# binary and sources
bin_PROGRAMS = gcritic2
gcritic2_SOURCES = main.cpp camera.cpp view.cpp shapes.cpp

# subdirectories
SUBDIRS=imgui .

# static libraries
gcritic2_LDADD=$(top_srcdir)/src/libcritic2.a \
	$(top_srcdir)/src/oldlibs/libmisc.a \
	$(top_srcdir)/src/ciftbx/libciftbx.a \
	$(top_srcdir)/src/cubpack/libcubpack.a \
	$(top_srcdir)/src/spglib/libspglib.a \
	$(top_srcdir)/src/qhull/libqhull_critic.a \
	imgui/imgui.o imgui/imgui_draw.o imgui/imgui_impl_glfw_gl3.o \
	imgui/imgui_demo.o imgui/imgui_dock.o imgui/imgui_widgets.o \
	imgui/gl3w.o

gcritic2_DEPENDENCIES=$(gcritic2_LDADD)

# Include the header files
AM_CXXFLAGS=-I$(top_srcdir)/include `pkg-config --cflags glfw3` -Wall -Wformat -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11 -w -g
AM_LDFLAGS=-ldl -lgfortran -lgomp -lglfw -lGLU -lGL
