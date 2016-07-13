#!/usr/bin/gnuplot -persist

set terminal x11 size 600,600
# set xrange [0:1]
# set yrange [0:1]
# set zrange [0:1]
# set xlabel "x-coordinate"
# set ylabel "y-coordinate"
# set zlabel "z-coordinate"
# set cblabel "Density"
# unset grid
# unset key

splot 'data/frame0.dat' u 1:2:3 with points pointtype 7 ps 1.4
