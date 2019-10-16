load "stylefile.gp"
set terminal qt 0
plot filename using "t":"T+R" title "Normalised (T+R)(E)"
set terminal qt 1
set xrange [0:8500]
set xlabel "t (fm/c)"
set log y
set format y '10^{%L}'
plot filename using "t":"N_0" w l title "Channel 1"
