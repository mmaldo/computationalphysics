set xlabel "x_n-1"
set ylabel "x_n"
set style data linespoints
set size square
set xrange [0:3]
unset xzeroaxis
set yrange [0:3]
unset yzeroaxis
plot "fixedPointsolve.data" using 2:3 title "Fixed Point Solve"

set term push  
set term postscript eps enh color solid lw 4 font 18
set out "fixedPointsolve.eps"
replot
set out 
set term pop
