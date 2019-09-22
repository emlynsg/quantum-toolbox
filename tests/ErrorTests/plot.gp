load "stylefile.gp"
set output "wavepacket.svg"
set xlabel 'x (a.u.)'
set key top right
set yrange [-0.2:0.5]

plot "<paste cN_x_dt_1.000000dx_1.000000.txt cN_norm_dt_1.000000dx_1.000000.txt" using ($1):($2) with lines title 'dt 0.01 dx 0.01'

