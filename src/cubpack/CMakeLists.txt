set(SOURCES buckley.f90 check.f90 cui.f90 divide.f90 ds_routines.f90
            error_handling.f90 global_all.f90 internal_types.f90 region_processor.f90
            rule_1.f90 rule_c2.f90 rule_c3.f90 rule_cn.f90 rule_general.f90
            rule_t2.f90 rule_t3.f90 rule_tn.f90 volume.f90)

add_library(cubpack OBJECT ${SOURCES})
