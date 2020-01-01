 set terminal postscript eps color enhanced size 10,7 'Helvetica' 40
 set style line 1 lt 1 lc rgbcolor variable lw 4
 set style line 2 lt 1 lc rgbcolor "#007700"
 set style line 3 lt 4 lc rgbcolor "#0000ff"
 set output "012_plane_mols_plane.eps"
 set size ratio -1
 unset key
 set xlabel "x (ang_ )"
 set ylabel "y (ang_ )"
 set xrange [  0.000:  6.000] 
 set yrange [  0.000:  6.000]
 load "012_plane_mols_plane-label.gnu"
plot "012_plane_mols_plane.iso" with lines ls 2
