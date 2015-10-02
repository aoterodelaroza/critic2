#! /bin/bash

## Convert from a cif file to quantum espresso input
## Use: cif2scf.sh file.cif file.scf.in

critic2 <<EOF
crystal $1
write ${1%cif}scf.in
EOF
