#!/usr/bin/gnuplot -persist

set terminal x11 size 1000,1000
set xrange [0:1]
set yrange [0:1]
set xlabel "x"
set ylabel "y"
# set legend "Fluid density"
unset grid
unset key
do for [i=1:500] {
  	plot 'data/frame'.i.'.dat' u 1:2:3 with points palette pointtype 7 ps 1.2
	set title 'Frame '.i.''
	pause 0.05E0
}
