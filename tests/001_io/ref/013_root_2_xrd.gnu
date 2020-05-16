set terminal postscript eps color enhanced "Helvetica" 25
set output "013_root_2_xrd.eps"

set label 1 "111" at 36.9434413,15.6314903 center rotate by 90 font "Helvetica,12"
set label 2 "200" at 42.9199702,102.0000000 center rotate by 90 font "Helvetica,12"
set label 3 "220" at 62.3149386,52.2193348 center rotate by 90 font "Helvetica,12"
set label 4 "311" at 74.7019408,9.7186193 center rotate by 90 font "Helvetica,12"
set label 5 "222" at 78.6432820,16.5864044 center rotate by 90 font "Helvetica,12"

set xlabel "2{/Symbol Q} (degrees)"
set ylabel "Intensity (arb. units)"
set xrange [5.0000000:90.0000000]
set style data lines
set grid
unset key
plot "013_root_2_xrd.dat" w lines

