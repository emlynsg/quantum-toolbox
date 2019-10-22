load "stylefile.gp"
set xrange [80:150]
set yrange [0.000000000000001:1]
set log y
set xlabel "E/V0"
set ylabel "T"
plot "N_1.csv" using "E":(1.0-$2) w l title "None", "N_2.csv" using "E":(1.0-$2) w l title "2", "N_3.csv" u "E":(1.0-$2) w l title "3", "N_4.csv" u "E":(1.0-$2) w l title "4", "N_10.csv" u "E":(1.0-$2) w l title "10" lw 3 lc black
