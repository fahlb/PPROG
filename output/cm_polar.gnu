set terminal pdf
set output "polar(cm)_plot.pdf"

set ylabel "d_{cm} [10^3 km]"

# plot Janus and Saturn
set title "Movement of Saturn and Janus"
plot "cm_data.dat" u ($1/(365.25*24*60*60)):($3/1e6) w l title "Janus" ,\
	 "cm_data.dat" u ($1/(365.25*24*60*60)):($4/1e6) w l title "Epimetheus"

# plot difference between distances to cm as function of time 

set title " difference of distances to cm of Janus and Epimetheus"

plot "cm_data.dat" u ($1/(365.25*24*60*60)):($5/1e3) title "difference"
