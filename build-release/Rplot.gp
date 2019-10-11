load "stylefile.gp"
set xrange [0.1:4.5]
set yrange [0:1]
plot "N_3.csv" using "E/V0":"R0/R" w l title "Channel 1", "N_3.csv" using "E/V0":"R1/R" w l title "Channel 2", "N_3.csv" using "E/V0":"R2/R" w l title "Channel 3"
