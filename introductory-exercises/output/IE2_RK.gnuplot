set terminal pdf size 6,6
set output "IE2_RK.pdf"

# plot x(t)
set title "position x as function of time t"
set xlabel "t [s]"
set ylabel "x(t) [m]"

plot "IE2_RK_h=0.100.dat" u 1:2 title "h=0.1s	" w l,\
	 "IE2_RK_h=0.010.dat" u 1:2 title "h=0.01s	" w l,\
	 "IE2_RK_h=0.001.dat" u 1:2 title "h=0.001s " w l,\
	 sin(x)

# plot v(t)
set title "velocity v as function of time t"
set xlabel "t [s]"
set ylabel "v(t) [m/s]"

plot "IE2_RK_h=0.100.dat" u 1:3 title "h=0.1s	" w l,\
	 "IE2_RK_h=0.010.dat" u 1:3 title "h=0.01s	" w l,\
	 "IE2_RK_h=0.001.dat" u 1:3 title "h=0.001s " w l,\
	 cos(x)

# plot x(t) against v(t)
set title "velocity as function of position"
set xlabel "x(t) [m]"
set ylabel "v(x(t)) [m/s]"
set size square

set object 1 circle at 0,0 size 1

plot "IE2_RK_h=0.100.dat" u 2:3 title "h=0.1s  " w l,\
	 "IE2_RK_h=0.010.dat" u 2:3 title "h=0.01s " w l,\
	 "IE2_RK_h=0.001.dat" u 2:3 title "h=0.001s" w l

