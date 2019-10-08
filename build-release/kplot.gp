set multiplot layout 2,2 rowsfirst
# --- GRAPH a
set label 1 'init' at graph 0.92,0.9 font ',12'
plot "DassoFig80.csv" using "k":"init" with lines ls 1
# --- GRAPH b
set label 1 '0' at graph 0.92,0.9 font ',12'
plot "DassoFig80.csv" using "k":"g" with lines, "DassoFig80.csv" using "k":"ex" with lines
# --- GRAPH c
set label 1 '2' at graph 0.92,0.9 font ',12'
plot "DassoFig82.csv" using "k":"g" with lines, "DassoFig82.csv" using "k":"ex" with lines
# --- GRAPH d
set label 1 '-2' at graph 0.92,0.9 font ',12'
plot "DassoFig8-2.csv" using "k":"g" with lines, "DassoFig8-2.csv" using "k":"ex" with lines
unset multiplot
