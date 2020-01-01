set encoding iso_8859_1
set terminal postscript eps color enhanced "Helvetica"
set output "010_plane_colormap_01-colormap.eps"

# line styles
set style line 1 lt 1 lw 1 lc rgb "#000000"
set style line 2 lt 1 lw 1 lc rgb "#000000"

# title, key, size
unset title
unset key
set size ratio -1

# set pm3d at b map interpolate 5,5
set pm3d at b map

# tics
set cbtics

# color schemes
set palette defined ( 0 "red", 1 "white", 2 "green" ) 

# set contours
unset clabel
set contour base
set cntrparam bspline
# set cntrparam levels incremental -min,step,max

load '010_plane_colormap_01-label.gnu'
splot "010_plane_colormap_01" u 4:5:(log(abs($6))) ls 1 w pm3d notitle
