set term pdfcairo
set outp "frankHertz.pdf"

set xlabel "[U] = V"
set ylabel "Anzahl an ankommenden Elektronen"

plot "output.dat" using 1:2 w l title "Simulation"

unset term
unset outp