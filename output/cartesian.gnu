set terminal pdf
set output "cartesian_plot.pdf"

set size square
set xlabel "x [10^3 km]"
set ylabel "y [10^3 km]"
set xrange [-200:200]
set yrange [-200:200]

# plot Janus and Saturn
set title "Movement of Saturn and Janus"
plot "data.dat" u ($2/1e6):($3/1e6) lc 3 title "Saturn",\
	 "data.dat" u ($6/1e6):($7/1e6) lc 1 title "Janus"

# plot Epimetheus and Saturn
set title "Movement of Saturn and Epimetheus"
plot "data.dat" u ($2/1e6):($3/1e6) lc 3 title "Saturn",\
	 "data.dat" u ($10/1e6):($11/1e6) lc 5 title "Epimetheus"

# plot Janus, Epimetheus and Saturn
set title "Movement of Saturn, Epimetheus and Janus"
plot "data.dat" u ($2/1e6):($3/1e6) lc 3 title "Saturn",\
	 "data.dat" u ($6/1e6):($7/1e6) lc 1 title "Janus" ,\
	 "data.dat" u ($10/1e6):($11/1e6) lc 5 title "Epimetheus"

