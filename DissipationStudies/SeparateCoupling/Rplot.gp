load "stylefile.gp"
set xrange [0.1:4.5]
set yrange [0:1]
set log x
set xlabel "E/V0"
set ylabel "R"
#plot "N_1.csv" using "E/V0":"N_R0" w l title "Channel 1"
#plot "N_2.csv" using "E/V0":"N_R0" w l title "Channel 1", "N_2.csv" #using "E/V0":"N_R1" w l title "Channel 2"
#plot "N_3.csv" using "E/V0":"N_R0" w l title "Channel 1", "N_3.csv" #using "E/V0":"N_R1" w l title "Channel 2", "N_3.csv" using "E/#V0":"N_R2" w l title "Channel 3"
plot "N_4.csv" using "E/V0":"N_R0" w l title "Channel 1", "N_4.csv" using "E/V0":"N_R1" w l title "Channel 2", "N_4.csv" using "E/V0":"N_R2" w l title "Channel 3", "N_4.csv" using "E/V0":"N_R3" w l title "Channel 4"
