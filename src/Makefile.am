# Critic2 targets and sources
dist_bin_SCRIPTS = 
bin_PROGRAMS = critic2
critic2_SOURCES = critic2.F90 cfftnd.f90 \
	spglib.f90 keast.f90 param.F90 spgs.f90 tools_system.c \
	sympg.f90 \
	doqhull.c \
	quadpack.f90\
	tricks.f90 \
	abinit_private.f90 abinit_private@proc.f90 \
	autocp.f90 autocp@proc.f90 \
	arithmetic.F90 arithmetic@proc.F90 \
	bader.f90 bader@proc.f90 \
	bisect.f90 bisect@proc.f90 \
	config.f90 config@proc.F90 \
	crystalseedmod.f90 crystalseedmod@proc.f90 \
	crystalmod.f90 crystalmod@proc.f90 \
	dftb_private.f90 dftb_private@proc.f90 \
	elk_private.f90 elk_private@proc.f90 \
	fieldmod.f90 fieldmod@proc.f90 \
	fieldseedmod.f90 fieldseedmod@proc.f90 \
	flux.f90 flux@proc.f90 \
	fragmentmod.f90 fragmentmod@proc.f90 \
	global.f90 global@proc.f90 \
	graphics.f90 graphics@proc.f90 \
	grid1mod.f90 grid1mod@proc.f90 \
	grid3mod.f90 grid3mod@proc.f90 \
	hashmod.f90 hashmod@proc.f90 \
	hirshfeld.f90 hirshfeld@proc.f90 \
	integration.f90 integration@proc.f90 \
	meshmod.f90 meshmod@proc.f90 \
	molcalc.f90 molcalc@proc.F90 \
	nci.f90 nci@proc.f90 \
	pi_private.f90 pi_private@proc.f90 \
	qtree_basic.f90 qtree_basic@proc.f90 \
	qtree.f90 qtree@proc.f90 \
	qtree_gpaths.f90 qtree_gpaths@proc.f90 \
	qtree_tetrawork.f90 qtree_tetrawork@proc.f90 \
	qtree_utils.f90 qtree_utils@proc.f90 \
	rhoplot.f90 rhoplot@proc.f90 \
	stm.f90 stm@proc.f90 \
	struct_drivers.f90 struct_drivers@proc.f90 \
	surface.f90 surface@proc.f90 \
	systemmod.f90 systemmod@proc.f90 \
	tools.f90 tools@proc.f90 \
	tools_io.f90 tools_io@proc.F90 \
	tools_math.f90 tools_math@proc.f90 tools_math@lebedev.f90 \
	types.f90 types@proc.f90 \
	wfn_private.f90 wfn_private@proc.f90 \
	wien_private.f90 wien_private@proc.f90 \
	xdm.f90 xdm@proc.f90 \
	yt.f90 yt@proc.f90

if GUI
  noinst_LIBRARIES=libcritic2.a
  libcritic2_a_SOURCES = gui_interface.f90 c_interface_module.f90 cfftnd.f90 \
	spglib.f90 keast.f90 param.F90 spgs.f90 tools_system.c \
	sympg.f90 \
	doqhull.c quadpack.f90\
	tricks.f90 \
	abinit_private.f90 abinit_private@proc.f90 \
	autocp.f90 autocp@proc.f90 \
	arithmetic.F90 arithmetic@proc.F90 \
	bisect.f90 bisect@proc.f90 \
	bader.f90 bader@proc.f90 \
	config.f90 config@proc.F90 \
	crystalseedmod.f90 crystalseedmod@proc.f90 \
	crystalmod.f90 crystalmod@proc.f90 \
	dftb_private.f90 dftb_private@proc.f90 \
	elk_private.f90 elk_private@proc.f90 \
	fieldmod.f90 fieldmod@proc.f90 \
	fieldseedmod.f90 fieldseedmod@proc.f90 \
	flux.f90 flux@proc.f90 \
	fragmentmod.f90 fragmentmod@proc.f90 \
	global.f90 global@proc.f90 \
	graphics.f90 graphics@proc.f90 \
	grid1mod.f90 grid1mod@proc.f90 \
	grid3mod.f90 grid3mod@proc.f90 \
	hashmod.f90 hashmod@proc.f90 \
	hirshfeld.f90 hirshfeld@proc.f90 \
	integration.f90 integration@proc.f90 \
	meshmod.f90 meshmod@proc.f90 \
	molcalc.f90 molcalc@proc.F90 \
	nci.f90 nci@proc.f90 \
	pi_private.f90 pi_private@proc.f90 \
	qtree_basic.f90 qtree_basic@proc.f90 \
	qtree.f90 qtree@proc.f90 \
	qtree_gpaths.f90 qtree_gpaths@proc.f90 \
	qtree_tetrawork.f90 qtree_tetrawork@proc.f90 \
	qtree_utils.f90 qtree_utils@proc.f90 \
	rhoplot.f90 rhoplot@proc.f90 \
	stm.f90 stm@proc.f90 \
	struct_drivers.f90 struct_drivers@proc.f90 \
	surface.f90 surface@proc.f90 \
	systemmod.f90 systemmod@proc.f90 \
	tools.f90 tools@proc.f90 \
	tools_io.f90 tools_io@proc.F90 \
	tools_math.f90 tools_math@proc.f90 tools_math@lebedev.f90 \
	types.f90 types@proc.f90 \
	wfn_private.f90 wfn_private@proc.f90 \
	wien_private.f90 wien_private@proc.f90 \
	xdm.f90 xdm@proc.f90 \
	yt.f90 yt@proc.f90

endif

# Library subdirectories
SUBDIRS = oldlibs ciftbx cubpack qhull spglib

# Static libraries
critic2_LDADD=oldlibs/libmisc.a ciftbx/libciftbx.a cubpack/libcubpack.a \
	spglib/libspglib.a $(LIBADD)
critic2_DEPENDENCIES=$(critic2_LDADD)

# Include the subdirectories with the modules
AM_FCFLAGS= -I$(top_srcdir)/src/cubpack

# Data directory
if !WIN
  AM_FCFLAGS+=-DDATADIR='"$(datadir)"'
else
  AM_FCFLAGS+=-DWIN
endif

# Compilation flags (in config.F90)
if !WIN
  DATE2 := $(shell date || echo "(config) @DATE@")
  AM_FCFLAGS+=-DDATE='"$(DATE2)"'
  COMMIT := $(shell git rev-parse --short HEAD || echo "(config) @VERSION@")
  AM_FCFLAGS+=-DCOMMIT='"$(COMMIT)"'
endif

# Clean the modules as well
mostlyclean-local:
	-rm -rf *.mod *.smod

# Modules depend on objects
%.o: %.f90
	@rm -f $@ $(@:.o=.mod) $(@:.o=.smod)
	$(AM_V_PPFC)$(PPFCCOMPILE) -c -o $@ $<
%.o: %.F90
	@rm -f $@ $(@:.o=.mod) $(@:.o=.smod)
	$(AM_V_PPFC)$(PPFCCOMPILE) -c -o $@ $<
%.mod: %.f90
	@rm -f $@ $(@:.mod=.o) $(@:.mod=.smod)
	$(AM_V_PPFC)$(PPFCCOMPILE) -c -o $(@:.mod=.o) $<
%.mod: %.F90
	@rm -f $@ $(@:.mod=.o) $(@:.mod=.smod)
	$(AM_V_PPFC)$(PPFCCOMPILE) -c -o $(@:.mod=.o) $<
%.smod: %.f90
	@rm -f $@ $(@:.smod=.o) $(@:.smod=.mod)
	$(AM_V_PPFC)$(PPFCCOMPILE) -c -o $(@:.smod=.o) $<
%.smod: %.F90
	@rm -f $@ $(@:.smod=.o) $(@:.smod=.mod)
	$(AM_V_PPFC)$(PPFCCOMPILE) -c -o $(@:.smod=.o) $<

# Always recompile the config file to get the compilation details right
config@proc.o : .FORCE

.FORCE: 

#xx# Dependencies: generated by makemake.sh
grid3mod@proc.o grid3mod@proc.smod crystalseedmod@proc.o crystalseedmod@proc.smod : abinit_private.mod
critic2.o critic2.mod : autocp.mod
integration@proc.o integration@proc.smod : bader.mod
critic2.o critic2.mod qtree@proc.o qtree@proc.smod : bisect.mod
gui_interface.o gui_interface.mod : c_interface_module.mod
critic2.o critic2.mod gui_interface.o gui_interface.mod global@proc.o global@proc.smod : config.mod
yt.o yt.mod xdm@proc.o xdm@proc.smod systemmod.o systemmod.mod systemmod.smod struct_drivers@proc.o struct_drivers@proc.smod struct_drivers.o struct_drivers.mod struct_drivers.smod meshmod@proc.o meshmod@proc.smod meshmod.o meshmod.mod meshmod.smod integration@proc.o integration@proc.smod fieldmod.o fieldmod.mod fieldmod.smod : crystalmod.mod
systemmod@proc.o systemmod@proc.smod systemmod.o systemmod.mod systemmod.smod struct_drivers@proc.o struct_drivers@proc.smod struct_drivers.o struct_drivers.mod struct_drivers.smod integration@proc.o integration@proc.smod gui_interface.o gui_interface.mod crystalmod@proc.o crystalmod@proc.smod crystalmod.o crystalmod.mod crystalmod.smod autocp@proc.o autocp@proc.smod : crystalseedmod.mod
fieldmod.o fieldmod.mod fieldmod.smod : dftb_private.mod
fieldmod.o fieldmod.mod fieldmod.smod : elk_private.mod
molcalc@proc.o molcalc@proc.smod critic2.o critic2.mod xdm@proc.o xdm@proc.smod systemmod@proc.o systemmod@proc.smod systemmod.o systemmod.mod systemmod.smod stm@proc.o stm@proc.smod rhoplot@proc.o rhoplot@proc.smod qtree@proc.o qtree@proc.smod nci@proc.o nci@proc.smod meshmod@proc.o meshmod@proc.smod meshmod.o meshmod.mod meshmod.smod integration@proc.o integration@proc.smod integration.o integration.mod integration.smod hirshfeld@proc.o hirshfeld@proc.smod bisect@proc.o bisect@proc.smod autocp@proc.o autocp@proc.smod : fieldmod.mod
systemmod@proc.o systemmod@proc.smod fieldmod.o fieldmod.mod fieldmod.smod : fieldseedmod.mod
critic2.o critic2.mod : flux.mod
nci@proc.o nci@proc.smod nci.o nci.mod nci.smod fieldmod@proc.o fieldmod@proc.smod fieldmod.o fieldmod.mod fieldmod.smod crystalmod@proc.o crystalmod@proc.smod crystalmod.o crystalmod.mod crystalmod.smod : fragmentmod.mod
critic2.o critic2.mod xdm@proc.o xdm@proc.smod systemmod@proc.o systemmod@proc.smod surface@proc.o surface@proc.smod struct_drivers@proc.o struct_drivers@proc.smod stm@proc.o stm@proc.smod rhoplot@proc.o rhoplot@proc.smod quadpack.o quadpack.mod qtree_utils@proc.o qtree_utils@proc.smod qtree_tetrawork@proc.o qtree_tetrawork@proc.smod qtree_gpaths@proc.o qtree_gpaths@proc.smod qtree_basic@proc.o qtree_basic@proc.smod qtree_basic.o qtree_basic.mod qtree_basic.smod qtree@proc.o qtree@proc.smod qtree.o qtree.mod qtree.smod pi_private@proc.o pi_private@proc.smod nci@proc.o nci@proc.smod integration@proc.o integration@proc.smod gui_interface.o gui_interface.mod grid1mod@proc.o grid1mod@proc.smod flux@proc.o flux@proc.smod fieldseedmod@proc.o fieldseedmod@proc.smod fieldmod@proc.o fieldmod@proc.smod crystalseedmod@proc.o crystalseedmod@proc.smod crystalmod@proc.o crystalmod@proc.smod bisect@proc.o bisect@proc.smod autocp@proc.o autocp@proc.smod : global.mod
integration@proc.o integration@proc.smod flux@proc.o flux@proc.smod crystalmod@proc.o crystalmod@proc.smod crystalmod.o crystalmod.mod crystalmod.smod bisect@proc.o bisect@proc.smod autocp@proc.o autocp@proc.smod : graphics.mod
critic2.o critic2.mod xdm@proc.o xdm@proc.smod struct_drivers@proc.o struct_drivers@proc.smod pi_private.o pi_private.mod pi_private.smod nci@proc.o nci@proc.smod hirshfeld@proc.o hirshfeld@proc.smod gui_interface.o gui_interface.mod fieldmod@proc.o fieldmod@proc.smod dftb_private.o dftb_private.mod dftb_private.smod crystalmod@proc.o crystalmod@proc.smod : grid1mod.mod
xdm@proc.o xdm@proc.smod rhoplot@proc.o rhoplot@proc.smod integration@proc.o integration@proc.smod hirshfeld@proc.o hirshfeld@proc.smod fieldmod@proc.o fieldmod@proc.smod fieldmod.o fieldmod.mod fieldmod.smod crystalmod@proc.o crystalmod@proc.smod crystalmod.o crystalmod.mod crystalmod.smod : grid3mod.mod
param.o param.mod arithmetic@proc.o arithmetic@proc.smod arithmetic.o arithmetic.mod arithmetic.smod systemmod.o systemmod.mod systemmod.smod grid3mod.o grid3mod.mod grid3mod.smod fieldmod@proc.o fieldmod@proc.smod fieldmod.o fieldmod.mod fieldmod.smod : hashmod.mod
critic2.o critic2.mod : hirshfeld.mod
critic2.o critic2.mod qtree@proc.o qtree@proc.smod bisect@proc.o bisect@proc.smod : integration.mod
qtree_utils@proc.o qtree_utils@proc.smod qtree@proc.o qtree@proc.smod : keast.mod
molcalc@proc.o molcalc@proc.smod xdm@proc.o xdm@proc.smod : meshmod.mod
critic2.o critic2.mod : molcalc.mod
critic2.o critic2.mod : nci.mod
fieldmod.o fieldmod.mod fieldmod.smod : pi_private.mod
critic2.o critic2.mod : qtree.mod
qtree_utils@proc.o qtree_utils@proc.smod qtree_utils.o qtree_utils.mod qtree_utils.smod qtree_tetrawork@proc.o qtree_tetrawork@proc.smod qtree_tetrawork.o qtree_tetrawork.mod qtree_tetrawork.smod qtree_gpaths@proc.o qtree_gpaths@proc.smod qtree_gpaths.o qtree_gpaths.mod qtree_gpaths.smod qtree@proc.o qtree@proc.smod : qtree_basic.mod
qtree_tetrawork@proc.o qtree_tetrawork@proc.smod : qtree_gpaths.mod
qtree@proc.o qtree@proc.smod : qtree_tetrawork.mod
qtree@proc.o qtree@proc.smod : qtree_utils.mod
integration@proc.o integration@proc.smod : quadpack.mod
critic2.o critic2.mod : rhoplot.mod
crystalmod@proc.o crystalmod@proc.smod crystalmod.o crystalmod.mod crystalmod.smod : spglib.mod
critic2.o critic2.mod gui_interface.o gui_interface.mod crystalseedmod@proc.o crystalseedmod@proc.smod : spgs.mod
critic2.o critic2.mod : stm.mod
critic2.o critic2.mod nci@proc.o nci@proc.smod autocp@proc.o autocp@proc.smod : struct_drivers.mod
qtree_basic@proc.o qtree_basic@proc.smod qtree@proc.o qtree@proc.smod integration@proc.o integration@proc.smod integration.o integration.mod integration.smod graphics@proc.o graphics@proc.smod graphics.o graphics.mod graphics.smod bisect@proc.o bisect@proc.smod bisect.o bisect.mod bisect.smod autocp@proc.o autocp@proc.smod : surface.mod
crystalmod@proc.o crystalmod@proc.smod : sympg.mod
molcalc@proc.o molcalc@proc.smod critic2.o critic2.mod yt.o yt.mod xdm@proc.o xdm@proc.smod struct_drivers@proc.o struct_drivers@proc.smod struct_drivers.o struct_drivers.mod struct_drivers.smod stm@proc.o stm@proc.smod rhoplot@proc.o rhoplot@proc.smod quadpack.o quadpack.mod qtree_utils@proc.o qtree_utils@proc.smod qtree_tetrawork@proc.o qtree_tetrawork@proc.smod qtree_gpaths@proc.o qtree_gpaths@proc.smod qtree_basic@proc.o qtree_basic@proc.smod qtree@proc.o qtree@proc.smod nci@proc.o nci@proc.smod integration@proc.o integration@proc.smod integration.o integration.mod integration.smod hirshfeld@proc.o hirshfeld@proc.smod gui_interface.o gui_interface.mod flux@proc.o flux@proc.smod bisect@proc.o bisect@proc.smod bisect.o bisect.mod bisect.smod bader@proc.o bader@proc.smod bader.o bader.mod bader.smod autocp@proc.o autocp@proc.smod : systemmod.mod
yt.o yt.mod integration@proc.o integration@proc.smod fieldmod@proc.o fieldmod@proc.smod crystalmod@proc.o crystalmod@proc.smod bisect@proc.o bisect@proc.smod autocp@proc.o autocp@proc.smod : tools.mod
molcalc@proc.o molcalc@proc.smod critic2.o critic2.mod arithmetic@proc.o arithmetic@proc.smod yt.o yt.mod xdm@proc.o xdm@proc.smod wien_private@proc.o wien_private@proc.smod wfn_private@proc.o wfn_private@proc.smod wfn_private.o wfn_private.mod wfn_private.smod tools_math@proc.o tools_math@proc.smod tools_math@lebedev.o tools_math@lebedev.smod tools_math.o tools_math.mod tools_math.smod tools@proc.o tools@proc.smod systemmod@proc.o systemmod@proc.smod systemmod.o systemmod.mod systemmod.smod sympg.o sympg.mod surface@proc.o surface@proc.smod struct_drivers@proc.o struct_drivers@proc.smod stm@proc.o stm@proc.smod spgs.o spgs.mod rhoplot@proc.o rhoplot@proc.smod qtree_utils@proc.o qtree_utils@proc.smod qtree_tetrawork@proc.o qtree_tetrawork@proc.smod qtree_gpaths@proc.o qtree_gpaths@proc.smod qtree_basic@proc.o qtree_basic@proc.smod qtree@proc.o qtree@proc.smod pi_private@proc.o pi_private@proc.smod nci@proc.o nci@proc.smod meshmod@proc.o meshmod@proc.smod integration@proc.o integration@proc.smod hirshfeld@proc.o hirshfeld@proc.smod gui_interface.o gui_interface.mod grid3mod@proc.o grid3mod@proc.smod grid1mod@proc.o grid1mod@proc.smod graphics@proc.o graphics@proc.smod global@proc.o global@proc.smod fragmentmod@proc.o fragmentmod@proc.smod flux@proc.o flux@proc.smod fieldseedmod@proc.o fieldseedmod@proc.smod fieldmod@proc.o fieldmod@proc.smod fieldmod.o fieldmod.mod fieldmod.smod elk_private@proc.o elk_private@proc.smod dftb_private@proc.o dftb_private@proc.smod crystalseedmod@proc.o crystalseedmod@proc.smod crystalmod@proc.o crystalmod@proc.smod bisect@proc.o bisect@proc.smod bader@proc.o bader@proc.smod bader.o bader.mod bader.smod autocp@proc.o autocp@proc.smod abinit_private@proc.o abinit_private@proc.smod : tools_io.mod
yt.o yt.mod xdm@proc.o xdm@proc.smod wien_private@proc.o wien_private@proc.smod systemmod@proc.o systemmod@proc.smod sympg.o sympg.mod surface@proc.o surface@proc.smod struct_drivers@proc.o struct_drivers@proc.smod rhoplot@proc.o rhoplot@proc.smod qtree_tetrawork@proc.o qtree_tetrawork@proc.smod qtree_basic@proc.o qtree_basic@proc.smod qtree@proc.o qtree@proc.smod pi_private@proc.o pi_private@proc.smod nci@proc.o nci@proc.smod meshmod@proc.o meshmod@proc.smod integration@proc.o integration@proc.smod grid3mod@proc.o grid3mod@proc.smod graphics@proc.o graphics@proc.smod fragmentmod@proc.o fragmentmod@proc.smod flux@proc.o flux@proc.smod fieldmod@proc.o fieldmod@proc.smod elk_private@proc.o elk_private@proc.smod dftb_private@proc.o dftb_private@proc.smod crystalseedmod@proc.o crystalseedmod@proc.smod crystalmod@proc.o crystalmod@proc.smod bisect@proc.o bisect@proc.smod bader@proc.o bader@proc.smod autocp@proc.o autocp@proc.smod : tools_math.mod
critic2.o critic2.mod : tricks.mod
molcalc@proc.o molcalc@proc.smod arithmetic@proc.o arithmetic@proc.smod arithmetic.o arithmetic.mod arithmetic.smod yt.o yt.mod xdm@proc.o xdm@proc.smod wfn_private@proc.o wfn_private@proc.smod wfn_private.o wfn_private.mod wfn_private.smod systemmod@proc.o systemmod@proc.smod systemmod.o systemmod.mod systemmod.smod sympg.o sympg.mod struct_drivers@proc.o struct_drivers@proc.smod rhoplot@proc.o rhoplot@proc.smod rhoplot.o rhoplot.mod rhoplot.smod qtree_tetrawork@proc.o qtree_tetrawork@proc.smod qtree_gpaths@proc.o qtree_gpaths@proc.smod qtree_gpaths.o qtree_gpaths.mod qtree_gpaths.smod qtree_basic@proc.o qtree_basic@proc.smod pi_private@proc.o pi_private@proc.smod pi_private.o pi_private.mod pi_private.smod nci@proc.o nci@proc.smod meshmod@proc.o meshmod@proc.smod integration@proc.o integration@proc.smod grid3mod@proc.o grid3mod@proc.smod grid3mod.o grid3mod.mod grid3mod.smod grid1mod@proc.o grid1mod@proc.smod graphics@proc.o graphics@proc.smod fragmentmod@proc.o fragmentmod@proc.smod fragmentmod.o fragmentmod.mod fragmentmod.smod flux@proc.o flux@proc.smod fieldseedmod@proc.o fieldseedmod@proc.smod fieldmod@proc.o fieldmod@proc.smod fieldmod.o fieldmod.mod fieldmod.smod dftb_private@proc.o dftb_private@proc.smod dftb_private.o dftb_private.mod dftb_private.smod crystalseedmod@proc.o crystalseedmod@proc.smod crystalseedmod.o crystalseedmod.mod crystalseedmod.smod crystalmod@proc.o crystalmod@proc.smod crystalmod.o crystalmod.mod crystalmod.smod bisect@proc.o bisect@proc.smod bader@proc.o bader@proc.smod bader.o bader.mod bader.smod autocp@proc.o autocp@proc.smod : types.mod
integration@proc.o integration@proc.smod fieldmod@proc.o fieldmod@proc.smod fieldmod.o fieldmod.mod fieldmod.smod crystalseedmod@proc.o crystalseedmod@proc.smod : wfn_private.mod
fieldmod.o fieldmod.mod fieldmod.smod : wien_private.mod
critic2.o critic2.mod : xdm.mod
integration@proc.o integration@proc.smod : yt.mod
critic2.o critic2.mod yt.o yt.mod systemmod@proc.o systemmod@proc.smod rhoplot@proc.o rhoplot@proc.smod grid3mod@proc.o grid3mod@proc.smod grid3mod.o grid3mod.mod grid3mod.smod global@proc.o global@proc.smod fieldmod@proc.o fieldmod@proc.smod fieldmod.o fieldmod.mod fieldmod.smod crystalseedmod@proc.o crystalseedmod@proc.smod bader@proc.o bader@proc.smod : arithmetic.mod
tools_io@proc.o tools_io@proc.smod molcalc@proc.o molcalc@proc.smod critic2.o critic2.mod arithmetic@proc.o arithmetic@proc.smod yt.o yt.mod xdm@proc.o xdm@proc.smod wien_private@proc.o wien_private@proc.smod wfn_private@proc.o wfn_private@proc.smod tools_math@proc.o tools_math@proc.smod tools_math@lebedev.o tools_math@lebedev.smod tools_io.o tools_io.mod tools_io.smod systemmod@proc.o systemmod@proc.smod sympg.o sympg.mod surface@proc.o surface@proc.smod struct_drivers@proc.o struct_drivers@proc.smod stm@proc.o stm@proc.smod spgs.o spgs.mod rhoplot@proc.o rhoplot@proc.smod quadpack.o quadpack.mod qtree_utils@proc.o qtree_utils@proc.smod qtree_basic@proc.o qtree_basic@proc.smod pi_private@proc.o pi_private@proc.smod nci@proc.o nci@proc.smod meshmod@proc.o meshmod@proc.smod integration@proc.o integration@proc.smod gui_interface.o gui_interface.mod grid3mod@proc.o grid3mod@proc.smod grid1mod@proc.o grid1mod@proc.smod graphics@proc.o graphics@proc.smod global@proc.o global@proc.smod global.o global.mod global.smod fragmentmod@proc.o fragmentmod@proc.smod flux@proc.o flux@proc.smod fieldseedmod@proc.o fieldseedmod@proc.smod fieldseedmod.o fieldseedmod.mod fieldseedmod.smod fieldmod@proc.o fieldmod@proc.smod fieldmod.o fieldmod.mod fieldmod.smod elk_private@proc.o elk_private@proc.smod dftb_private@proc.o dftb_private@proc.smod crystalseedmod@proc.o crystalseedmod@proc.smod crystalmod@proc.o crystalmod@proc.smod crystalmod.o crystalmod.mod crystalmod.smod config.o config.mod config.smod bisect@proc.o bisect@proc.smod bader@proc.o bader@proc.smod autocp@proc.o autocp@proc.smod : param.mod
