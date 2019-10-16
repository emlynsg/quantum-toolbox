load "stylefile.gp"
set terminal pdfcairo size 20cm,15cm
set xrange [0:700]
set xlabel "t (fm/c)"
set log y
set format y '10^{%L}'
set key bottom right
Ns = "1 2 3 4"
Es = "0.5 0.9 1 1.2 2"

do for [i=1:words(Ns)] {
do for [j=1:words(Es)] {
filename = 'N_'.word(Ns, i).'_E_'.word(Es,j).'.csv'
set output 'N_'.word(Ns, i).'_E_'.word(Es,j).'.pdf'
plot for [k=1:word(Ns, i)] filename using "t":"N_".(k-1) w l title "Channel ".(k-1)
}
}
