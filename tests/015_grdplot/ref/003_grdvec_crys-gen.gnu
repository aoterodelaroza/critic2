 set terminal postscript eps color enhanced size 10,7 'Helvetica' 40
 set style line 1 lt 1 lc rgbcolor variable lw 4
 set style line 2 lt 1 lc rgbcolor "#007700"
 set style line 3 lt 4 lc rgbcolor "#0000ff"
 set output "003_grdvec_crys-gen.eps"
 set size ratio -1
 unset key
 set xlabel "x (bohr )"
 set ylabel "y (bohr )"
 set xrange [  0.000:  8.753] 
 set yrange [  0.000:  8.753]
 load "003_grdvec_crys-gen-label.gnu"
plot "003_grdvec_crys-gen.iso" with lines ls 2, "003_grdvec_crys-gen.dat" u 1:2:6 with lines ls 1
