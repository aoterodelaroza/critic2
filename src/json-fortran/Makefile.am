noinst_LIBRARIES=libjsonfortran.a

libjsonfortran_a_SOURCES=json_file_module.F90 json_kinds.F90 json_module.F90 json_parameters.F90 json_string_utilities.F90 json_value_module.F90

json_module.o json_module.mod : json_file_module.mod
json_value_module.o json_value_module.mod json_string_utilities.o json_string_utilities.mod json_parameters.o json_parameters.mod json_module.o json_module.mod json_file_module.o json_file_module.mod : json_kinds.mod
json_value_module.o json_value_module.mod json_string_utilities.o json_string_utilities.mod json_module.o json_module.mod json_file_module.o json_file_module.mod : json_parameters.mod
json_value_module.o json_value_module.mod json_module.o json_module.mod json_file_module.o json_file_module.mod : json_string_utilities.mod
json_module.o json_module.mod json_file_module.o json_file_module.mod : json_value_module.mod

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

# clean the modules as well
mostlyclean-local:
	-rm -rf *.$(MODEXT)

