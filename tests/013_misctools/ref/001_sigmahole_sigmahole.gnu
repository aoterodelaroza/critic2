set terminal wxt size 800,800
set mapping cartesian

set style line 1 lt 1 lw 1 lc rgb "#000000"
set style line 2 lt 1 lw 1 lc rgb "#000000"
unset clabel

set palette defined (-0.03 "red", 0 "white", 0.03 "blue")

unset title
unset key

set view equal xyz
set xlabel "x (ang)"
set ylabel "y (ang)"
set zlabel "z (ang)"
set view 45
set pm3d depthorder # interpolate 0,0
set cbrange [-0.03:0.03]

splot "001_sigmahole_sigmahole.dat" u 5:6:7:9 with pm3d
pause -1
