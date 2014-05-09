set xlabel "n"
set ylabel "|x_n-sqrt(5)|"
set style data linespoints
set logscale y
set xrange[1:10]
set size square
plot "fixedPointsolve.data" using 1:4 title "Fixed Point Convergence"

set term push  
set term postscript eps enh color solid lw 4 font 18
set out "fixedPointconvergence.eps"
replot
set out 
set term pop
