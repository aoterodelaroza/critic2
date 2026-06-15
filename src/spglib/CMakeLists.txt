set(SOURCES arithmetic.c cell.c debug.c delaunay.c determination.c
  hall_symbol.c kgrid.c kpoint.c magnetic_spacegroup.c mathfunc.c
  msg_database.c niggli.c overlap.c pointgroup.c primitive.c refinement.c
  sitesym_database.c site_symmetry.c spacegroup.c spg_database.c spglib.c
  spglib_f.c spin.c symmetry.c
  )

add_library(spglib OBJECT ${SOURCES})

# spglib >= 2.x targets the C11 standard; the top-level project sets C99, so
# bump just this object library.
set_target_properties(spglib PROPERTIES C_STANDARD 11)
