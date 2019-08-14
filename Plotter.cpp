#include "Plotter.h"
#include "System.h"

Plotter::Plotter(System sys, const int &ncols, const bool &showpsi, const bool &showpotential,
                 const bool &showenergy, const bool &shownorm, const bool &showavgx, const bool &showpsik)
                 : system(Wavefunction(Grid(1, 0.0, 1.0, 1.0), 1.0), Potential(Grid(1, 0.0, 1.0, 1.0))) {
  system = sys;
  nColumns = ncols;
  showPsi = showpsi;
  showPotential = showpotential;
  showEnergy = showenergy;
  showNorm = shownorm;
  showAvgX = showavgx;
  showPsiK = showpsik;

}

Plotter::~Plotter(){
  std::cout << "Plotter deleted" << std::endl;
}

void Plotter::test(){
  std::cout << "Test Plotter" << std::endl;
  GnuplotPipe gp;
  setPlotStyle(gp, 0);
  gp.sendLine(R"(set output "test.pdf")");
  gp.sendLine(R"(set xlabel "x axis label")");
  gp.sendLine(R"(set ylabel "y axis label")");
  gp.sendLine(R"(set key top right)");
  gp.sendLine(R"(plot [-pi/2:pi] cos(x),-(sin(x) > sin(x+1) ? sin(x) : sin(x+1)))");
  gp.sendLine(R"(set terminal pop)");
  gp.sendLine(R"(set output)");
  gp.sendLine(R"(replot)");
}

void Plotter::plot(){

}

void Plotter::setPlotStyle(GnuplotPipe &g, int stylenum){
  if (stylenum == 1) {
    g.sendLine(R"(cd "..")");
    g.sendLine(R"(set encoding utf8)");
    g.sendLine(R"(set style line 1 lc rgb '#E41A1C' pt 1 ps 1 lt 1 lw 2)");
    g.sendLine(R"(set style line 2 lc rgb '#377EB8' pt 6 ps 1 lt 1 lw 2)");
    g.sendLine(R"(set style line 3 lc rgb '#4DAF4A' pt 2 ps 1 lt 1 lw 2)");
    g.sendLine(R"(set style line 4 lc rgb '#984EA3' pt 3 ps 1 lt 1 lw 2)");
    g.sendLine(R"(set style line 5 lc rgb '#FF7F00' pt 4 ps 1 lt 1 lw 2)");
    g.sendLine(R"(set style line 6 lc rgb '#FFFF33' pt 5 ps 1 lt 1 lw 2)");
    g.sendLine(R"(set style line 7 lc rgb '#A65628' pt 7 ps 1 lt 1 lw 2)");
    g.sendLine(R"(set style line 8 lc rgb '#F781BF' pt 8 ps 1 lt 1 lw 2)");
    g.sendLine(R"(set palette maxcolors 8)");
    g.sendLine(R"(set palette defined ( 0 '#E41A1C', 1 '#377EB8', 2 '#4DAF4A', 3 '#984EA3',\
4 '#FF7F00', 5 '#FFFF33', 6 '#A65628', 7 '#F781BF' ))");
    g.sendLine(R"(set style line 11 lc rgb '#808080' lt 1 lw 3)");
    g.sendLine(R"(set border 0 back ls 11)");
    g.sendLine(R"(set tics out nomirror)");
    g.sendLine(R"(set style line 12 lc rgb '#808080' lt 0 lw 1)");
    g.sendLine(R"(set grid back ls 12)");
    g.sendLine(R"(set terminal pdfcairo enhanced color dashed font "Alegreya, 14" \
rounded size 16 cm, 9.6 cm)");
  }
  else {
    g.sendLine(R"(cd "..")");
    g.sendLine(R"(set terminal pdfcairo font "Gill Sans,7" linewidth 3 rounded fontscale 1.0)");
    g.sendLine(R"(set style line 80 lt rgb "#808080")");
    g.sendLine(R"(set style line 81 lt 0)");
    g.sendLine(R"(set style line 81 lt rgb "#808080")");
    g.sendLine(R"(set grid back linestyle 81)");
    g.sendLine(R"(set border 3 back linestyle 80)");
    g.sendLine(R"(set xtics nomirror)");
    g.sendLine(R"(set ytics nomirror)");
    g.sendLine(R"(set style line 1 lt rgb "#A00000" lw 2 pt 1)");
    g.sendLine(R"(set style line 2 lt rgb "#00A000" lw 2 pt 6)");
    g.sendLine(R"(set style line 3 lt rgb "#5060D0" lw 2 pt 2)");
    g.sendLine(R"(set style line 4 lt rgb "#F25900" lw 2 pt 9)");
  }

}