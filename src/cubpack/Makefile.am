noinst_LIBRARIES=libcubpack.a

libcubpack_a_SOURCES=buckley.f90 check.f90 cui.f90 divide.f90 ds_routines.f90 error_handling.f90 global_all.f90 internal_types.f90 region_processor.f90 rule_1.f90 rule_c2.f90 rule_c3.f90 rule_cn.f90 rule_general.f90 rule_t2.f90 rule_t3.f90 rule_tn.f90 volume.f90

# dependencies
ds_routines.o: buckley.o
rule_1.o: buckley.o
rule_c2.o: buckley.o
rule_t2.o: buckley.o
rule_c3.o: buckley.o
rule_t3.o: buckley.o
rule_tn.o: buckley.o
rule_cn.o: buckley.o
volume.o: buckley.o
divide.o: buckley.o internal_types.o
rule_general.o: rule_t3.o rule_t2.o rule_c2.o rule_c3.o rule_1.o \
	        internal_types.o
region_processor.o : divide.o rule_general.o
global_all.o: region_processor.o volume.o ds_routines.o internal_types.o
cui.o: global_all.o check.o error_handling.o ds_routines.o internal_types.o
internal_types.o: buckley.o
check.o: internal_types.o
rule_general.o: rule_tn.o rule_cn.o

# clean the modules as well
mostlyclean-local:
	-rm -rf *.$(MODEXT)

