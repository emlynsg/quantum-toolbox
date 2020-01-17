load "stylefile.gp"
set output "transmission.svg"
set xlabel 'E/V0'
set key bottom right
set yrange [0.0:1.2]
set xrange[0.0:6.0]
plot "<paste SQ_Barrier_Erel.txt SQ_Barrier_T_V0_100.txt" using ($1):($2) with lines title 'V0=100MeV'

#a = 1.0
#hbar = 197.3
#mu = 1.0
#T(x) = 
