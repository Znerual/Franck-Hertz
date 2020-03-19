set term pdfcairo
set outp "frankHertz.pdf"

set title "Franck-Hertz Versuch"

set xlabel "[U] = V"
set ylabel "[I] = mA"


plot "output.dat" using 1:($2*1000) w l title ""

unset term
unset outp