set xrange [3.3:3.6]
plot "results.data" using 2:1
set term pop
set terminal png
set out "logmap3.png"
replot
