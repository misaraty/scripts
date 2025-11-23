set terminal pdfcairo color enhanced font "Times,25" fontscale 0.4 size 4,3 lw 1
set ylabel "{/=30 {Max. eigenvalue}" offset 1, 0
set xlabel "{/=30 {Temperature (K)}}" offset 0, 0.5
set lmargin screen 0.18
set rmargin screen 0.96
set tmargin screen 0.94
set bmargin screen 0.20
set tics font ",30"
set xtics 1
set ytics 2
set mxtics 2
set mytics 2
set key font ",20"
set key spacing 1.5
set xrange [0:6]
set yrange [0:6]
set out "fig4.pdf"
set key at graph 0.9, 0.9
set arrow  from 0, 1 to 6, 1 nohead
set arrow  from 4.61, 2.5 to 4.61, 1.2
set label 1 "T_c = 4.61 K" font ",25"
set label 1 at graph 0.65, 0.5
plot "data_max_eigenvalue.dat" with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle
reset
