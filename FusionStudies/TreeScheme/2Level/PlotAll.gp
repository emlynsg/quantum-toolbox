load "stylefile.gp"
#set terminal pdfcairo size 20cm,15cm
#set xrange [0:5000]
set xlabel "t (fm/c)"
#set log y
#set format y '10^{%L}'
set key bottom right
#Ns = "2 3 4 5 6 7 8"

#do for [i=1:words(Ns)] {
#filename = "N_".word(Ns, i)."ChannelsOverTime.csv"
#set output "N_".word(Ns, i)."ChannelsOverTime.pdf"
#plot for [k=1:word(Ns, i)] filename using "t":"Norm_".(k-1) w l title #"Channel ".(k-1)
#}

#plot "N_2ChannelsOverTime.csv" using "t":"Norm_0" w l title "Level 1", #"N_2ChannelsOverTime.csv" using "t":"Norm_1" w l title "Level 2 #Channel 1","N_2ChannelsOverTime.csv" using "t":"Norm_2" w l title #"Level 2 Channel 2"

plot "N_2ChannelsOverTime_Tree.csv" using "t":"Norm_L0" w l title "Level 1", "N_2ChannelsOverTime_Tree.csv" using "t":"Norm_L1" w l title "Level 2"
