set terminal pdfcairo color enhanced font "Times,25" fontscale 0.4 size 4,3 lw 1
set ylabel "{/=30 {Δ_{nk} (meV)}" offset 1, 0
set xlabel "{/=30 {Temperature (K)}}" offset 0, 0.5
set lmargin screen 0.18
set rmargin screen 0.96
set tmargin screen 0.94
set bmargin screen 0.20
set tics font ",30"
set xtics 10
set ytics 3
set mxtics 5
set mytics 1
set key font ",18"
set key samplen 2
set key at graph 1.0, 0.90
set xrange [0:60]
set yrange [0:15]
set out "fig11.pdf"
scale=3E-3
set label 1 "μ_c^*=0.05" font ",18"
set label 1 at graph 0.80, 0.95
set style fill transparent solid 0.2
plot \
for [temp in "10 15 20 25 30 35 40 45 50 55"] gprintf('../epw1-FSR/mgb2.imag_aniso_gap0_%06.2f',temp) u (($3)+scale*($5)):($2) with filledc above x1=temp lc rgb "red" lw 0 notitle, \
for [temp in "10"] gprintf('../epw1-FSR/mgb2.imag_aniso_gap0_%06.2f',temp) u (($3)+scale*($5)):($2) w l lc rgb "red" dt 3 lw 2 title "FSR", \
for [temp in "15 20 25 30 35 40 45 50 55"] gprintf('../epw1-FSR/mgb2.imag_aniso_gap0_%06.2f',temp) u (($3)+scale*($5)):($2) w l lc rgb "red" dt 3 lw 2 notitle, \
for [temp in "10 15 20 25 30 35 40 45 50 55"] gprintf('./mgb2.imag_aniso_gap0_%06.2f',temp) u (($3)+scale*($5)):($2) with filledc above x1=temp lc rgb "blue" lw 0 notitle, \
for [temp in "10"] gprintf('./mgb2.imag_aniso_gap0_%06.2f',temp) u (($3)+scale*($5)):($2) w l lc rgb "blue" lw 1 title "FBW", \
for [temp in "15 20 25 30 35 40 45 50 55"] gprintf('./mgb2.imag_aniso_gap0_%06.2f',temp) u (($3)+scale*($5)):($2) w l lc rgb "blue" lw 1 notitle
reset