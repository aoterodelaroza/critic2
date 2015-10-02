#! /bin/bash

## charges.sh: 

CRITIC2="critic2"
usage='\n
   charges.sh: calculate atomic charges in a density grid using critic2\n
\n
   Usage:\n
\t    charges.sh file_DEN\n
\t    charges.sh file.cube\n
\t    charges.sh CHGCAR (reads CONTCAR if it exists, or POSCAR otherwise. Also, the POTCAR)\n
\t    charges.sh -d ... # dry run (print the input file but do not run it)\n
\t    charges.sh -n ... # search for non-nuclear maxima
\t    charges.sh -h # this message\n
\n
'    

if [ ! $(which $CRITIC2) ] ; then
    echo "critic2 executable not found"
    exit 1
fi

if [ $# != 1 ] ; then
    echo -e $usage
    exit 1
fi

tempfile="temp."$$".incritic"
outfile="temp."$$".scf.in"
if [ ${1##*_} == "DEN" ] || [ ${1##*.} == "cube" ]; then
    cat > $tempfile <<EOF
noguess
crystal $1
write $outfile noprimitive
EOF
elif [ ${1} == "POSCAR" ] || [ ${1} == "CONTCAR" ] ; then
    # the atoms
    if [ ! -f "POTCAR" ] ; then
	echo "POTCAR not found"
	exit 1
    fi	
    # the input file
    cat > $tempfile <<EOF
noguess
crystal $1 POTCAR
write $outfile noprimitive
EOF
fi
$CRITIC2 $tempfile >& /dev/null
rm -f $tempfile

