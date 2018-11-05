#! /bin/bash

mpirun -np 4 pw.x < feo_paw.scf.in | tee feo_paw.scf.out
mpirun -np 4 pp.x < feo.rhoae.in
mpirun -np 4 pw.x < feo.scf.in | tee feo.scf.out
mpirun -np 4 pp.x < feo.rho.in

mpirun -np 4 open_grid.x < feo.opengrid.in | tee feo.opengrid.out
pw2critic.x < feo.pw2critic.in

cat > feo_up.win <<EOF
num_wann = 9
num_bands = 9
num_iter = 20000
conv_tol = 1e-6
conv_window = 3

begin unit_cell_cart
bohr
   -0.000000000000     4.093146803830    -4.093146803830
    4.093146803830     4.093146803830     0.000000000000
    4.093146803830     0.000000000000    -4.093146803830
end unit_cell_cart

begin atoms_frac
fe    0.00000000     0.00000000     0.00000000
o     0.50000000     0.50000000     0.50000000
end atoms_frac

begin projections
random
end projections

search_shells = 24
kmesh_tol = 1d-3
mp_grid : 6 6 6

begin kpoints
EOF
awk '/List to be put/,/^ *$/' feo.opengrid.out | grep -v List | grep -v '^ *$' >> feo_up.win
echo "end kpoints" >> feo_up.win

mpirun -np 4 wannier90.x -pp feo_up.win
mpirun -np 4 pw2wannier90.x < feo.pw2wan.up.in | tee feo.pw2wan.up.out
mpirun -np 4 wannier90.x feo_up.win

cat > feo_dn.win <<EOF
num_wann = 5
num_bands = 9
exclude_bands : 6-9
num_iter = 20000
conv_tol = 1e-6
conv_window = 3

begin unit_cell_cart
bohr
   -0.000000000000     4.093146803830    -4.093146803830
    4.093146803830     4.093146803830     0.000000000000
    4.093146803830     0.000000000000    -4.093146803830
end unit_cell_cart

begin atoms_frac
fe    0.00000000     0.00000000     0.00000000
o     0.50000000     0.50000000     0.50000000
end atoms_frac

begin projections
random
end projections

search_shells = 24
kmesh_tol = 1d-3
mp_grid : 6 6 6

begin kpoints
EOF
awk '/List to be put/,/^ *$/' feo.opengrid.out | grep -v List | grep -v '^ *$' >> feo_dn.win
echo "end kpoints" >> feo_dn.win

mpirun -np 4 wannier90.x -pp feo_dn.win
mpirun -np 4 pw2wannier90.x < feo.pw2wan.dn.in | tee feo.pw2wan.dn.out

cat > feo_dn.win <<EOF
num_wann = 5
num_bands = 5
num_iter = 20000
conv_tol = 1e-6
conv_window = 3

begin unit_cell_cart
bohr
   -0.000000000000     4.093146803830    -4.093146803830
    4.093146803830     4.093146803830     0.000000000000
    4.093146803830     0.000000000000    -4.093146803830
end unit_cell_cart

begin atoms_frac
fe    0.00000000     0.00000000     0.00000000
o     0.50000000     0.50000000     0.50000000
end atoms_frac

begin projections
random
end projections

search_shells = 24
kmesh_tol = 1d-3
mp_grid : 6 6 6

begin kpoints
EOF
awk '/List to be put/,/^ *$/' feo.opengrid.out | grep -v List | grep -v '^ *$' >> feo_dn.win
echo "end kpoints" >> feo_dn.win

mpirun -np 4 wannier90.x feo_dn.win
