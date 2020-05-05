set terminal postscript eps color enhanced "Helvetica" 25
set output "016_rdf_options_gen.eps"

set xlabel "r (bohr)"
set ylabel "RDF(r)"
set xrange [1.0000000:23.0000000]
set style data lines
set grid
unset key
plot "016_rdf_options_gen.dat" w lines

