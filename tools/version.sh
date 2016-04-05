#! /bin/sh

command -v git >/dev/null 2>&1 || { echo "no-git"; exit; }

if [ -d .git ] ; then
    git rev-parse --short HEAD
elif [ -f ChangeLog ] ; then
    head -n 1 ChangeLog | cut -f2 -d' '
else
    echo "unknown"
fi
