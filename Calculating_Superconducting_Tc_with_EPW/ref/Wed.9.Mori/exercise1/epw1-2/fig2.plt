set terminal pdfcairo color enhanced font "Times,25" fontscale 0.4 size 4,3 lw 1
set ylabel "{/=30 {Δ_0 (meV)}" offset 1, 0
set xlabel "{/=30 {Temperature (K)}}" offset 0, 0.5
set lmargin screen 0.18
set rmargin screen 0.96
set tmargin screen 0.94
set bmargin screen 0.20
set tics font ",30"
set xtics 1
set ytics 0.4
set mxtics 2
set mytics 4
set key font ",20"
set xrange [0:6]
set yrange [0:1.2]
set out "fig2.pdf"
set key at graph 0.9, 0.9
plot "pb.imag_iso_gap0" with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle
reset