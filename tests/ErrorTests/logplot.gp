load "stylefile.gp"
# set output "wavepacket_log.svg"
set xlabel 'x (a.u.)'
set key top right
set yrange [*:10.0]
set log y
set format y '10^{%L}'

plot "<paste cN_x_dt_1.000000dx_1.000000.txt cN_norm_dt_1.000000dx_1.000000.txt" using ($1):($2) with lines title 'dt 0.01 dx 0.01'
