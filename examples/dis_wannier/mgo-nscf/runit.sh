#! /bin/bash

mpirun -np 4 pw.x < mgo_paw.scf.in | tee mgo_paw.scf.out
mpirun -np 4 pp.x < mgo.rhoae.in
mpirun -np 4 pw.x < mgo.scf.in | tee mgo.scf.out
mpirun -np 4 pp.x < mgo.rho.in
mpirun -np 4 pw.x < mgo.nscf.in | tee mgo.nscf.out

pw2critic.x < mgo.pw2critic.in

wannier90.x -pp mgo.win
mpirun -np 4 pw2wannier90.x < mgo.pw2wan.in
wannier90.x mgo.win

