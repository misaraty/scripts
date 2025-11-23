set terminal pdfcairo color enhanced font "Times,25" fontscale 0.4 size 7,6 lw 1
set out "fig5.pdf"
set multiplot
set border lw 3
dl1=0.1
dl2=0.1
dl3=0.101015254455
dl4=0.096225044865
dl5=0.096423651980 
length1=1.0
length2=1.5
length3=2.207106781187
length4=3.073132184971
length5=4.133792356751

set ylabel "{/=30 {f_{nest}(q)    (arb. units)}" offset 1, 0
set lmargin screen 0.07
set rmargin screen 0.48
set tmargin screen 0.94
set bmargin screen 0.55
set tics font ",30"
set ytics 2
set mytics 2
set key font ",20"
set key spacing 1.5
set xrange [0:length5]
set yrange [0:5]
set grid xtics lw 2 lt 1 lc 0
set xtics ("Γ"  0, "X"  length1, "W"  length2, "L"  length3, "Γ"  length4, "K"  length5)
unset key
plot \
"pb.nesting_fn" every ::0::10   u (($1-1) *dl1+0):($2/100)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"pb.nesting_fn" every ::10::15 u (($1-11)*dl2+length1):($2/100)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"pb.nesting_fn" every ::15::22 u (($1-16)*dl3+length2):($2/100)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"pb.nesting_fn" every ::22::31 u (($1-23)*dl4+length3):($2/100)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"pb.nesting_fn" every ::31::42 u (($1-32)*dl5+length4):($2/100)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle#, \
#"pb.nesting_fn" every ::42::42 u (($1-43)*0  +length5):($2/100)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle

set ylabel "{/=30 {λ_{qν} (ν=1)}" offset 1, 0
set lmargin screen 0.57
set rmargin screen 0.98
set tmargin screen 0.94
set bmargin screen 0.55
set tics font ",30"
set ytics 1
set mytics 2
set key font ",20"
set key spacing 1.5
set xrange [0:length5]
set yrange [0:2.8]
set grid xtics lw 2 lt 1 lc 0
set xtics ("Γ"  0, "X"  length1, "W"  length2, "L"  length3, "Γ"  length4, "K"  length5)
unset key
plot \
"lambda.phself.300.000K" every ::0::10   u (($1-1) *dl1+0):($2)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::10::15 u (($1-11)*dl2+length1):($2)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::15::22 u (($1-16)*dl3+length2):($2)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::22::31 u (($1-23)*dl4+length3):($2)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::31::42 u (($1-32)*dl5+length4):($2)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle#, \
#"lambda.phself.300.000K" every ::42::42 u (($1-43)*0  +length5):($2)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle

set ylabel "{/=30 {λ_{qν} (ν=2)}" offset 1, 0
set lmargin screen 0.07
set rmargin screen 0.48
set tmargin screen 0.44
set bmargin screen 0.05
set tics font ",30"
set ytics 1
set mytics 2
set key font ",20"
set key spacing 1.5
set xrange [0:length5]
set yrange [0:2.8]
set grid xtics lw 2 lt 1 lc 0
set xtics ("Γ"  0, "X"  length1, "W"  length2, "L"  length3, "Γ"  length4, "K"  length5)
unset key
plot \
"lambda.phself.300.000K" every ::0::10   u (($1-1) *dl1+0):($3)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::10::15 u (($1-11)*dl2+length1):($3)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::15::22 u (($1-16)*dl3+length2):($3)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::22::31 u (($1-23)*dl4+length3):($3)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::31::42 u (($1-32)*dl5+length4):($3)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle#, \
#"lambda.phself.300.000K" every ::42::42 u (($1-43)*0  +length5):($3)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle

set ylabel "{/=30 {λ_{qν} (ν=3)}" offset 1, 0
set lmargin screen 0.57
set rmargin screen 0.98
set tmargin screen 0.44
set bmargin screen 0.05
set tics font ",30"
set ytics 1
set mytics 2
set key font ",20"
set key spacing 1.5
set xrange [0:length5]
set yrange [0:2.8]
set grid xtics lw 2 lt 1 lc 0
set xtics ("Γ"  0, "X"  length1, "W"  length2, "L"  length3, "Γ"  length4, "K"  length5)
unset key
plot \
"lambda.phself.300.000K" every ::0::10   u (($1-1) *dl1+0):($4)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::10::15 u (($1-11)*dl2+length1):($4)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::15::22 u (($1-16)*dl3+length2):($4)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::22::31 u (($1-23)*dl4+length3):($4)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle, \
"lambda.phself.300.000K" every ::31::42 u (($1-32)*dl5+length4):($4)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle#, \
#"lambda.phself.300.000K" every ::42::42 u (($1-43)*0  +length5):($4)  with lp pt 7 ps 1 lw 2 dt 2 lc rgb "blue" notitle

reset