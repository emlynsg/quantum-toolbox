load "stylefile.gp"
# set output "wavepacket_log.svg"
set xlabel 'x (a.u.)'
set key top right
set yrange [*:10.0]
set log y
set format y '10^{%L}'

plot "<paste ErrorTest_x_dt_0.01_dx_0.01.txt ErrorTest_norm_dt_0.01_dx_0.01.txt" using ($1):($2) with lines title 'dt 0.01 dx 0.01', "<paste ErrorTest_x_dt_1_dx_1.txt ErrorTest_norm_dt_1_dx_1.txt" using ($1):($2) with lines title 'dt 1 dx 1', "<paste ErrorTest2_x_dt_1_dx_1.txt ErrorTest2_norm_dt_1_dx_1.txt" using ($1):($2) with lines title '2 dt 1 dx 1'
