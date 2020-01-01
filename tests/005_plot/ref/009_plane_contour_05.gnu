 set terminal postscript eps color enhanced size 10,7 'Helvetica' 40
 set style line 1 lt 1 lc rgbcolor variable lw 4
 set style line 2 lt 1 lc rgbcolor "#007700"
 set style line 3 lt 4 lc rgbcolor "#0000ff"
 set output "009_plane_contour_05.eps"
 set size ratio -1
 unset key
 set xlabel "x (bohr )"
 set ylabel "y (bohr )"
 set xrange [  0.000:  7.958] 
 set yrange [  0.000:  7.958]
 load "009_plane_contour_05-label.gnu"
plot "009_plane_contour_05.iso" with lines ls 2
