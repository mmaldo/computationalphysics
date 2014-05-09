set xrange [2.9:3.1]
plot "results.data" using 2:1
set term pop
set terminal png
set out "logmap2.png"
replot
