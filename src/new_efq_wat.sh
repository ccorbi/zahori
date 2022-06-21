#!/bin/sh
gnuplot <<EOF
set logscale x
#set xrange[0:$2]
set yrange[0:0.5]
set style line 1 lt 1 lw 2 pt 3


#p '2f0r_A_WAT$1.edf.dat' w l, '2f0r_B_WAT$1.edf.dat' w l,'3obs_WAT$1.edf.dat' w l,'3obq_WAT$1.edf.dat' w l,'3obu_WAT$1.edf.dat' w l,'3obx_WAT$1.edf.dat' w l,'4efe_ac_WAT$1.edf.dat' w l, '4eje_b_WAT$1.edf.dat' w l 

#p  '1s1q.ddbb_WAT$1' w l, '1s1q.cd.ddbb_WAT$1' w l,'2f0r_a.ddbb_WAT$1' w l,'2f0r_b.ddbb_WAT$1' w l,'4eje_AC.ddbb_WAT$1' w l, '4eje_BD.ddbb_WAT$1' w l

p  'random_0.01' w l, 'random_0.05' w l,'random_0.1' w l,'random_0.2' w l ,'random_0.3' w l, 'random_0.35' w l, 'random_0.4' w l, 'random_0.5' w l

set size square
set terminal postscript eps color
set output "plot_$1_$2_clusters_freq_log.ps"
replot

EOF
