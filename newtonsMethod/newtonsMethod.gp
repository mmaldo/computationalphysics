set xlabel "x_n-1"
set ylabel "x_n"
set xrange [0:3.5]
unset xzeroaxis
set yrange [0:3.5]
unset yzeroaxis
set style data linespoints
set size square
plot "newtonsMethod.data" using 2:3 title "Newton's Method"

set term push  
set term postscript eps enh color solid lw 4 font 18
set out "newtonsMethod.eps"
replot
set out 
set term pop
