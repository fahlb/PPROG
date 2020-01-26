set terminal pdf size 6,6
set output "IE1.pdf"


set xlabel "t [s]"
set ylabel "x(t) [m]"

plot "IE1_h=0.50.dat" u 1:2 title "h=0.5s ",\
	 "IE1_h=0.10.dat" u 1:2 title "h=0.1s ",\
	 "IE1_h=0.01.dat" u 1:2 title "h=0.01s",\
	 exp(-x)

