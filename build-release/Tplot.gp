load "stylefile.gp"
set xrange [50:150]
set yrange [0:1]
plot "N_2.csv" using "E":"T" w l title "2", "N_3.csv" u "E":"T" w l title "3", "N_4.csv" u "E":"T" w l title "4"
