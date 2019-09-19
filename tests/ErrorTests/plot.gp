load "stylefile.gp"
set output "wavepacket.svg"
set xlabel 'x (a.u.)'
set key top right
set yrange [-0.2:0.5]

plot "<paste ErrorTest_x_dt_0.01_dx_0.01.txt ErrorTest_norm_dt_0.01_dx_0.01.txt" using ($1):($2) with lines title '$\left\|\psi(x)\right\|$', "<paste ErrorTest_x_dt_0.01_dx_0.01.txt ErrorTest_real_dt_0.01_dx_0.01.txt" using ($1):($2) with lines title '$\Re(\psi(x))$', "<paste ErrorTest_x_dt_0.01_dx_0.01.txt ErrorTest_imag_dt_0.01_dx_0.01.txt" using ($1):($2) with lines title '$\Im(\psi(x))$'

