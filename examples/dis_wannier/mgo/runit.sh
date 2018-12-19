#! /bin/bash

mpirun -np 4 pw.x < mgo_paw.scf.in | tee mgo_paw.scf.out
mpirun -np 4 pp.x < mgo.rhoae.in
mpirun -np 4 pw.x < mgo.scf.in | tee mgo.scf.out
mpirun -np 4 pp.x < mgo.rho.in

mpirun -np 4 open_grid.x < mgo.opengrid.in | tee mgo.opengrid.out
pw2critic.x < mgo.pw2critic.in

cat > mgo.win <<EOF
num_wann = 4
num_iter = 20000
conv_tol = 1e-8
conv_window = 3

begin unit_cell_cart
bohr
  -0.000000000   3.901537427  -3.901537427
   3.901540775   3.901536617  -0.000000811
   3.901540775   0.000000811  -3.901536617
end unit_cell_cart

begin atoms_frac
mg       0.500000000   0.000000000   0.000000000
o        0.000000000   0.500000000   0.500000000
end atoms_frac

begin projections
random
end projections

search_shells = 24
kmesh_tol = 1d-3
mp_grid : 4 4 4

begin kpoints
EOF
awk '/List to be put/,/^ *$/' mgo.opengrid.out | grep -v List | grep -v '^ *$' >> mgo.win
echo "end kpoints" >> mgo.win

wannier90.x -pp mgo.win
mpirun -np 4 pw2wannier90.x < mgo.pw2wan.in
wannier90.x mgo.win


