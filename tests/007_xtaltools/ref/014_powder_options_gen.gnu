set terminal postscript eps color enhanced "Helvetica" 25
set output "014_powder_options_gen.eps"

set label 1 "001" at 18.4276715,102.0000000 center rotate by 90 font "Helvetica,12"

set xlabel "2{/Symbol Q} (degrees)"
set ylabel "Intensity (arb. units)"
set xrange [10.0000000:21.9500000]
set style data lines
set grid
unset key
plot "014_powder_options_gen.dat" w lines

