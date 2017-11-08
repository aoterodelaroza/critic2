noinst_PROGRAMS = imgui
imgui_SOURCES = imgui.cpp imgui_draw.cpp imgui_impl_glfw_gl3.cpp imgui_demo.cpp imgui_dock.cpp\
	imgui_widgets.cpp gl3w.c
imgui_LINK = true

AM_CXXFLAGS=-I$(top_srcdir)/include `pkg-config --cflags glfw3` -Wall -Wformat -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11 -w -g
