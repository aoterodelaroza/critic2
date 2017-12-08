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

dry=""
nnm=""
OPTSTRING=${@}
OPTIND=1
while getopts "dn" opt
do
	case $opt in
	    d)  dry="1" ;;
	    n)  nnm="nnm";;
	    ?)
		echo -e $usage
		exit 0
		;;
	esac
done
shift $((OPTIND-1))
if [ $# != 1 ] ; then
    echo -e $usage
    exit 1
fi

tempfile="temp."$$".incritic"
if [ ${1##*_} == "DEN" ] || [ ${1##*.} == "cube" ]; then
    cat > $tempfile <<EOF
crystal $1
load $1
yt $nnm
EOF
elif [ ${1} == "CHGCAR" ] || [ ${1} == "CHG" ] || [ ${1} == "AECCAR0" ] || [ ${1} == "AECCAR2" ] ; then
    # the structure
    if [ -f "CONTCAR" ] ; then
	struct="CONTCAR"
    elif [ -f "POSCAR" ] ; then
	struct="POSCAR"
    else
	echo "POSCAR/CONTCAR not found"
	exit 1
    fi
    # the atoms
    if [ ! -f "POTCAR" ] ; then
	echo "POTCAR not found"
	exit 1
    fi	
    atoms=$(grep VRHFIN POTCAR | cut -f2 -d= | cut -f1 -d: | sed 's/\\n/ /' | awk '{printf "%s ",$1}')
    # the input file
    cat > $tempfile <<EOF
crystal $struct $atoms
load $1
yt $nnm
EOF
fi
if [ -z $dry ]; then
    $CRITIC2 $tempfile
else
    cat $tempfile
fi
rm -f $tempfile

