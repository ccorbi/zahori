#!/bin/sh
gnuplot <<EOF
set xrange[0:$1]

p 'random_WATprob70_edf.dat' w l linecolor rgb 'black', 'random_WATprob80_edf.dat' w l linecolor rgb 'blue','random_WATprob90_edf.dat' w l linecolor rgb 'red','random_WATprob95_edf.dat' w l linecolor rgb 'green','random_WATprob99_edf.dat' w l linecolor rgb 'brown','random_WATprob995_edf.dat' w l linecolor rgb 'orange'

set size square
set terminal postscript eps color
set output "plot_random_ed_waters_scale_$1.ps"
replot

EOF
