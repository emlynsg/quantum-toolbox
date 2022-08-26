set key top right box opaque
set datafile separator ","
set key autotitle columnhead

# line styles
set linetype 1 dt 1 lc rgb 'black' lw 2#
set linetype 2 dt 3 lc rgb 'blue' lw 2#
set linetype 3 dt 4 lc rgb 'red' lw 2#
set linetype 4 dt "-" lc rgb '#41AB5D' lw 2#


# Standard border
set style line 11 lc rgb '#808080' lt 1 lw 3
set border 0 back ls 11
set tics out nomirror

# Standard grid
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12
unset grid

set grid
set xtics nomirror
set ytics nomirror
