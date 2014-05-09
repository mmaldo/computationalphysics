plot "results.data" using 2:1
set term push  
set term postscript eps enh color solid lw 4 font 18
set out "logmap.eps"
replot
set out 

