set xrange [0:10]
set xlabel "n"
set ylabel "Adjoint Eigenfunction"
plot for [j=2:11] "adj.out" u 1:j w l title 'Z_{'.j.'}'
