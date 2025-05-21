set style data dots
set nokey
set xrange [0: 6.07914]
set yrange [ -1.54120 : 19.64863]
set arrow from  1.28085,  -1.54120 to  1.28085,  19.64863 nohead
set arrow from  1.92128,  -1.54120 to  1.92128,  19.64863 nohead
set arrow from  2.82698,  -1.54120 to  2.82698,  19.64863 nohead
set arrow from  3.61134,  -1.54120 to  3.61134,  19.64863 nohead
set arrow from  4.96988,  -1.54120 to  4.96988,  19.64863 nohead
set xtics ("G"  0.00000,"X"  1.28085,"W"  1.92128,"L"  2.82698,"K"  3.61134,"G"  4.96988,"L"  6.07914)
 plot "pb_band.dat"
