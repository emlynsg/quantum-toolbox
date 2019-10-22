load "stylefile.gp"
set xrange [70:130]
set yrange [0.0000001:1]
set log y
set xlabel "E/V0"
set ylabel "T"
plot "N_1.csv" using "E":"T" w l title "None", "N_2.csv" using "E":"T" w l title "2", "N_3.csv" u "E":"T" w l title "3", "N_4.csv" u "E":"T" w l title "4", "N_10.csv" u "E":"T" w l title "10" lc black lt 3
