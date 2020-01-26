set terminal pdf size 6,6
set output "IE1.pdf"


set xlabel "t"
set ylabel "x(t)"

plot "IE1_h=0.10.dat", "IE1_h=0.50.dat", "IE1_h=0.01.dat", exp(-x) 1c 3

