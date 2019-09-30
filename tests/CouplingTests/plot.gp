load "stylefile.gp"
set xrange [4.00:16.00]
set yrange [0.0:1.5]
plot "DassoFig20.csv" using 1:6 with lines, "DassoFig22.csv" using 1:6 with lines
