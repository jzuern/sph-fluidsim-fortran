set xlabel "x"
set ylabel "y"
m="./data"
set terminal x11 0
set nokey
set grid
set title 'The parabola'
plot m using 1:2 with linespoints
