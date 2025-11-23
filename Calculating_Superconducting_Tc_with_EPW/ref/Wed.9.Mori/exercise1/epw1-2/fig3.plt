set terminal pdfcairo color enhanced font "Times,25" fontscale 0.4 size 4,3 lw 1
set ylabel "{/=30 {Δ_0 (meV)}" offset 1, 0
set xlabel "{/=30 {Temperature (K)}}" offset 0, 0.5
set lmargin screen 0.18
set rmargin screen 0.96
set tmargin screen 0.94
set bmargin screen 0.20
set tics font ",25"
set xtics 1
set ytics 0.4
set mxtics 2
set mytics 4
set key font ",18"
set xrange [0:6]
set yrange [0:1.2]
set out "fig3.pdf"
set key at graph 0.97, 0.95
plot  \
"pb.imag_iso_gap0" with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" title "Δ(iω_{min}) (lowest Matsubara freq.)", \
"pb.pade_iso_gap0" with lp pt 7 ps 1 lw 2 dt 2 lc rgb "black" title "Δ(ω=0) (Pade approx.)", \
"pb.acon_iso_gap0" with lp pt 7 ps 1 lw 2 dt 2 lc rgb "red" title "Δ(ω=0) (analytic cont.)"
reset