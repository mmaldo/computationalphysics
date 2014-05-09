set polar
set sample 2000
set zeroaxis
set trange [0:3.5*pi]
set xrange [-5:5]
set yrange [-5:5]
set size square
#set style line 1 lt 1 lc rgb "red" lw 3
#set style line 2 lt 1 lc rgb "blue" lw 3
#set style line 3 lt 1 lc rgb "black" lw 3
r(a,t) = exp(a*t)
x(k,t) = k*t
f(m,j,t)= sqrt(4*m*t)+ j
y(m,j,t)= -sqrt(4*m*t)+ j
a = .14
k = .29
m = .15
j = 0
plot r(a,t) title "Logarithmic Spiral" , \
x(k,t) title "Spiral of Arcimedes" , \
f(m,j,t)  title "Fermat Spiral" with lines 3, y(m,j,t) notitle with lines 3

set term push  
set term postscript eps enh color solid lw 4 font 18
set out "hw1p1.eps"
replot
set out 
set term pop