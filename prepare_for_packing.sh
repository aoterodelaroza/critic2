#! /bin/bash

make distclean
autoreconf --force -i
./configure
make
cd doc
./compile.sh
cd ..
