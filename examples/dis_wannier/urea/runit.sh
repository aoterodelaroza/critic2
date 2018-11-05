#! /bin/bash

mpirun -np 4 pw.x < urea_paw.scf.in | tee urea_paw.scf.out
mpirun -np 4 pp.x < urea.rhoae.in
mpirun -np 4 pw.x < urea.scf.in | tee urea.scf.out
mpirun -np 4 pp.x < urea.rho.in
mpirun -np 4 open_grid.x < urea.opengrid.in | tee urea.opengrid.out
pw2critic.x < urea.pw2critic.in

cat > urea.win <<EOF
num_wann = 24
num_iter = 20000
conv_tol = 1e-8
conv_window = 3

begin unit_cell_cart
bohr
   10.516325929510     0.000000000000     0.000000000000
    0.000000000000    10.516325929510     0.000000000000
    0.000000000000     0.000000000000     8.851477206440
end unit_cell_cart

begin atoms_frac
C    0.00000000     0.50000000     0.32600000
C    0.50000000     0.00000000     0.67400000
O    0.00000000     0.50000000     0.59530000
O    0.50000000     0.00000000     0.40470000
N    0.14590000     0.64590000     0.17660000
N    0.35410000     0.14590000     0.82340000
N    0.85410000     0.35410000     0.17660000
N    0.64590000     0.85410000     0.82340000
H    0.25750000     0.75750000     0.28270000
H    0.24250000     0.25750000     0.71730000
H    0.74250000     0.24250000     0.28270000
H    0.75750000     0.74250000     0.71730000
H    0.14410000     0.64410000     0.96200000
H    0.35590000     0.14410000     0.03800000
H    0.85590000     0.35590000     0.96200000
H    0.64410000     0.85590000     0.03800000
end atoms_frac

begin projections
random
end projections

search_shells = 24
kmesh_tol = 1d-3
mp_grid : 2 2 2

begin kpoints
EOF
awk '/List to be put/,/^ *$/' urea.opengrid.out | grep -v List | grep -v '^ *$' >> urea.win
echo "end kpoints" >> urea.win

mpirun -np 4 wannier90.x -pp urea.win
mpirun -np 4 pw2wannier90.x < urea.pw2wan.in
mpirun -np 4 wannier90.x urea.win

