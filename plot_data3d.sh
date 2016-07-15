#!/usr/bin/gnuplot -persist

set terminal x11 size 1000,1000
set xrange [0:1]
set yrange [0:1]
set zrange [0:1]
set xlabel "x-coordinate "
set ylabel "y-coordinate "
set zlabel "z-coordinate "
set border 4095
set grid

nfiles = system("find data/ -type f | wc -l")
do for [i=1:nfiles-1] {
  splot 'data/frame'.i.'.dat' u 1:2:3:4 with points palette pointtype 7 ps 2.0

	set title 'Frame '.i.''
	pause 0.2E0
}
