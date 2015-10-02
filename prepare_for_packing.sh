#! /bin/bash

rm ChangeLog
make ChangeLog
make distclean
autoreconf --force -i
./configure
make
cd doc
./compile.sh
cd ..
