#! /bin/sh

if [ -d .git ] ; then
    git rev-parse --short HEAD
elif [ -f ChangeLog ] ; then
    head -n 1 ChangeLog | cut -f2 -d' '
else
    echo "unknown version"
fi
