#!/usr/bin/gnuplot -persist

set terminal x11 size 1000,1000
set xrange [0:1]
set yrange [0:1]
set xlabel "x-coordinate"
set ylabel "y-coordinate"
set cblabel "Density"
unset grid
unset key

nfiles = system("find data/ -type f | wc -l")
do for [i=1:nfiles-1] {
  plot 'data/frame'.i.'.dat' u 1:2:3 with points palette pointtype 7 ps 1.4
	set title 'Frame '.i.''
	pause 0.05E0
}
