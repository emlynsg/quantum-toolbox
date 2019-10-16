load "stylefile.gp"
set xrange [70:130]
set yrange [0.000000001:1]
set xlabel "E/V0"
set ylabel "T"
plot "N_1.csv" using "E":"T" w l title "None", "N_2.csv" using "E":"T" w l title "2", "N_3.csv" u "E":"T" w l title "3", "N_4.csv" u "E":"T" w l title "4"
