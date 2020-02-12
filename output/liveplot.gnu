# use to track progress of running calculation
set terminal qt

set xrange [-0.4:1.2]
set yrange [0:13000]

plot "data.dat" u 2:3
pause 0.1
reread
