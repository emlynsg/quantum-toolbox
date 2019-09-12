load "stylefile.gp"
set output "transmission.svg"
set xlabel 'E/V0'
set key bottom right
set yrange [0.0:1.2]
plot "<paste SQ_Barrier_Energies_V0_100.txt SQ_Barrier_Transmissions_V0_100.txt" using ($1/100.0):($2) with lines title 'V0=100MeV', "<paste SQ_Barrier_Energies_V0_75.txt SQ_Barrier_Transmissions_V0_75.txt" using ($1/75.0):($2) with lines title 'V0=75MeV', "<paste SQ_Barrier_Energies_V0_40.txt SQ_Barrier_Transmissions_V0_40.txt" using ($1/40.0):($2) with lines title 'V0=40MeV', "<paste SQ_Barrier_Energies_V0_10.txt SQ_Barrier_Transmissions_V0_10.txt" using ($1/10.0):($2) with lines title 'V0=10MeV'

#a = 1.0
#hbar = 197.3
#mu = 1.0
#T(x) = 
