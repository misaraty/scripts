set terminal pdfcairo color enhanced font "Times,25" fontscale 0.4 size 4,3 lw 1
set ylabel "{/=30 {Δ (meV)}" offset 1, 0
set xlabel "{/=30 {Temperature (K)}}" offset 0, 0.5
set lmargin screen 0.18
set rmargin screen 0.99
set tmargin screen 0.99
set bmargin screen 0.20
set tics font ",30"
set xtics 2.0
set border lw 3
set key font ",14"
set xrange [0.0:16.6]
set yrange [0.0:3.8]
set out "fig15.pdf"
scale=8E-5
set ytics offset 0.5, 0
set key at graph 0.4, 0.15
set style fill transparent solid 0.2
plot \
for [temp in "0.2 6.0 10.0 12.0 13.0 13.4"] gprintf('./nb.imag_aniso_gap0_%06.2f',temp) u (($3)+scale*($5)):($2) with filledc above x1=temp lc rgb "blue" lt 1 lw 0 notitle, \
for [temp in "0.2"] gprintf('./nb.imag_aniso_gap0_%06.2f',temp) u (($3)+scale*($5)):($2) w l lc rgb "blue" lt 1 lw 1 title "μ_c^* = 0.25", \
for [temp in "6.0 10.0 12.0 13.0 13.4"] gprintf('./nb.imag_aniso_gap0_%06.2f',temp) u (($3)+scale*($5)):($2) w l lc rgb "blue" lt 1 lw 1 notitle
reset