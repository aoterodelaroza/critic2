set(SOURCES arithmetic.c cell.c debug.c delaunay.c determination.c
  hall_symbol.c kgrid.c kpoint.c mathfunc.c niggli.c overlap.c
  pointgroup.c primitive.c refinement.c sitesym_database.c
  site_symmetry.c spacegroup.c spg_database.c spglib.c spglib_f.c spin.c
  symmetry.c tetrahedron_method.c
  )

add_library(spglib OBJECT ${SOURCES})
