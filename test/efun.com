set xrange [0:10]
set xlabel "n"
set ylabel "Eigenfunction"
plot for [j=2:11] "efun.out" u 1:j w l title 'U_{'.j.'}'
