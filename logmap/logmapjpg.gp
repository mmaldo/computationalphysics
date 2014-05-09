set xlabel "r"
set ylabel "xn_1"
plot "results.data" using 2:1
set term pop
set terminal png
set out "logmap.png"
replot
