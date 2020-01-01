set terminal pdfcairo
set output "011_plane_relief_plane-relief.pdf"
set encoding iso_8859_1

set style line 1 lt 1 lc rgb "#000000" 

# Define the zrange and the capping functions
zmin = 0.00000E+00 
zmax = 2.00000E+00 
stats "011_plane_relief_plane.dat" u 6 nooutput
min(x) = (x<zmin)?min=x:zmin
max(x) = (x>zmax)?max=zmax:x
set zrange [(zmin<STATS_min)?STATS_min:zmin:(zmax>STATS_max)?STATS_max:zmax]

# tics, etc
unset colorbox
unset title
set format x "%.1f"
set format y "%.1f"

# Surface definition
set pm3d depthorder hidden3d 1
set hidden3d
set style fill transparent solid 0.7
set palette rgb 9,9,3
set view 60,45
set size ratio -1

splot "011_plane_relief_plane.dat" u 4:5:(max($6)) ls 1 w pm3d notitle
