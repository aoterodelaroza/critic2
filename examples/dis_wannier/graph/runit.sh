#! /bin/bash

mpirun -np 4 pw.x < graph_paw.scf.in | tee graph_paw.scf.out
mpirun -np 4 pp.x < graph.rhoae.in
mpirun -np 4 pw.x < graph.scf.in | tee graph.scf.out
mpirun -np 4 pp.x < graph.rho.in
mpirun -np 4 open_grid.x < graph.opengrid.in > graph.opengrid.out
pw2critic.x < graph.pw2critic.in > graph.pw2critic.out

cat > graph.win <<EOF
num_wann = 8
num_iter = 20000
conv_tol = 1e-8
conv_window = 3

begin unit_cell_cart
bohr
    4.641167382367     0.000000000000     0.000000000000
   -2.320583691184     4.019368856346     0.000000000000
    0.000000000000     0.000000000000    12.653606185802
end unit_cell_cart

begin atoms_frac
c    0.00000000     0.00000000     0.25000000
c    0.00000000     0.00000000     0.75000000
c    0.33333333     0.66666667     0.25000000
c    0.66666667     0.33333333     0.75000000
end atoms_frac

begin projections
random
end projections

search_shells = 24
kmesh_tol = 1d-3
mp_grid : 8 8 2

begin kpoints
EOF
awk '/List to be put/,/^ *$/' graph.opengrid.out | grep -v List | grep -v '^ *$' >> graph.win
echo "end kpoints" >> graph.win

mpirun -np 4 wannier90.x -pp graph.win
mpirun -np 4 pw2wannier90.x < graph.pw2wan.in > graph.pw2wan.out
mpirun -np 4 wannier90.x graph.win

