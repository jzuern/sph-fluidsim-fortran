#!/bin/bash
dir=data/*
for f in $dir; do
	gnuplot <<- EOF
		set terminal x11 size 1000,1000
		set xrange [0:1]
		set yrange [0:1]
		unset grid
		plot '${f}' u 1:2:3 with points palette pointtype 7 ps 1.2 title 'Particle visualization'
		pause 0.100E+00
EOF
done
