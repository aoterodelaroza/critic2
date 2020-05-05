set terminal postscript eps color enhanced "Helvetica" 25
set output "015_rdf_basic_rdf.eps"

set xlabel "r (bohr)"
set ylabel "RDF(r)"
set xrange [0.0000000:25.0000000]
set style data lines
set grid
unset key
plot "015_rdf_basic_rdf.dat" w lines

