set terminal pdfcairo color enhanced font "Times,25" fontscale 0.4 size 8,3 lw 1
set out "fig1.pdf"
set multiplot
set border lw 2
set key spacing 1.5


set lmargin screen 0.09
set rmargin screen 0.48
set tmargin screen 0.94
set bmargin screen 0.20
set ylabel "{/=30 {Δ (meV)}" offset 1, 0
set xlabel "{/=30 {iω (meV)}}" offset 0, 0.5
set tics font ",30"
set xtics 20
set ytics 0.4
set mxtics 2
set mytics 2
set key font ",20"
set xrange [0:60]
set yrange [-0.4:1.4]
plot "pb.imag_iso_000.30" u ($1*1000):($3*1000) w l lw 2 lt 1 lc rgb "black" notitle

reset
set border lw 2

set lmargin screen 0.59
set rmargin screen 0.98
set tmargin screen 0.94
set bmargin screen 0.20
set ylabel "{/=30 {Δ (meV)}" offset 1, 0
set xlabel "{/=30 {ω (meV)}}" offset 0, 0.5
set tics font ",30"
set xtics 20
set ytics 1.5
set mxtics 2
set mytics 3
set key font ",20"
set xrange [0:60]
set yrange [-2.5:3.5]
set key at graph 0.95, 0.9
plot "pb.pade_iso_000.30" u ($1*1000):($4*1000) w l lw 2 dt 1 lc rgb "black" title "Re(Δ)-Pade approx.", \
     "pb.pade_iso_000.30" u ($1*1000):($5*1000) w l lw 2 dt 2 lc rgb "red" title "Im(Δ)-Pade approx.", \
     "pb.acon_iso_000.30" u ($1*1000):($4*1000) w l lw 2 dt 4 lc rgb "green" title "Re(Δ)-analytic cont.", \
     "pb.acon_iso_000.30" u ($1*1000):($5*1000) w l lw 2 dt 3 lc rgb "blue" title "Im(Δ)-analytic cont."

unset multiplot
reset