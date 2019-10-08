load "stylefile.gp"
set xrange [90:105]
set log y
set yrange [0.0001:10]
plot "DassoFig80.csv" using "E":"T" with lines title "0", "DassoFig82.csv" using "E":"T" with lines title "2", "DassoFig8-2.csv" using "E":"T" with lines title "-2"
