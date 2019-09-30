load "stylefile.gp"
set xrange [0.37:2.0]
plot "DeltaT_0.01.csv" using 1:4 with lines, "DeltaT_0.01.csv" using 1:5 with lines
set term qt 1
plot "DeltaT_0.01.csv" using 1:($4-$5) with lines, "DeltaT_0.01.csv" using 1:3 with lines
set term qt 2
plot "DeltaT_0.01.csv" using 1:(($4-$5)/($5)) with lines
