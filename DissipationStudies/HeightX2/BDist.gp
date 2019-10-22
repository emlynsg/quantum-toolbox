load "stylefile.gp"
set xrange [70:200]
set log y
set xlabel "E/V0"
set ylabel "Barrier Distribution"
x0=NaN
y0=NaN
plot "N_1.csv" using (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w l title "None", "N_2.csv" using (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w l title "2", "N_3.csv" using (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w l title "3", "N_4.csv" using (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w l title "4"
