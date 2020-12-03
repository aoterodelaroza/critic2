set encoding iso_8859_1
set terminal postscript eps color enhanced 'Helvetica, 20'
set output '002_stm_current_stm.eps'

set palette rgb 34,35,36
set style line 1 lt 1 lw 1 lc rgb '#000000'
set style line 2 lt 1 lw 1 lc rgb '#000000'

set title 'Constant current plot, 5.000E-05 density

set xrange [0.000000E+00:1.156855E+01]
set yrange [0.000000E+00:1.644683E+01]
set xtics nomirror out
set ytics nomirror out
set xlabel 'x (\305)'
set ylabel 'y (\305)'

unset key
set size ratio -1

# set pm3d at b map interpolate 5,5
set pm3d at b map

set cbtics
set cblabel 'Distance to the surface (\305)'

splot '002_stm_current_stm.dat' u 3:4:5 w pm3d notitle ls 1
