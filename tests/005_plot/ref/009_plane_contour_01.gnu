 set terminal postscript eps color enhanced size 10,7 'Helvetica' 40
 set style line 1 lt 1 lc rgbcolor variable lw 4
 set style line 2 lt 1 lc rgbcolor "#007700"
 set style line 3 lt 4 lc rgbcolor "#0000ff"
 set output "009_plane_contour_01.eps"
 set size ratio -1
 unset key
 set xlabel "x (bohr )"
 set ylabel "y (bohr )"
 set xrange [  0.000: 11.549] 
 set yrange [  0.000: 11.549]
 load "009_plane_contour_01-label.gnu"
plot "009_plane_contour_01.iso" with lines ls 2, "009_plane_contour_01.neg.iso" with lines ls 3
