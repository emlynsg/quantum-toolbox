#include "Plotter.h"
#include "System.h"

Plotter::Plotter(System sys)
                 : system(Wavefunction(Grid(1, 0.0, 1.0, 1.0), 1.0), Potential(Grid(1, 0.0, 1.0, 1.0))) {
  system = sys;
}

Plotter::~Plotter(){
  std::cout << "Plotter deleted" << std::endl;
}

void Plotter::setPlotStyle(Gnuplot &g, int stylenum){
  if (stylenum == 1) {
    g << "cd '..'\n";
    g << "set encoding utf8\n";
    g << "set style line 1 lc rgb '#E41A1C' pt 1 ps 1 lt 1 lw 2\n";
    g << "set style line 2 lc rgb '#377EB8' pt 6 ps 1 lt 1 lw 2\n";
    g << "set style line 3 lc rgb '#4DAF4A' pt 2 ps 1 lt 1 lw 2\n";
    g << "set style line 4 lc rgb '#984EA3' pt 3 ps 1 lt 1 lw 2\n";
    g << "set style line 5 lc rgb '#FF7F00' pt 4 ps 1 lt 1 lw 2\n";
    g << "set style line 6 lc rgb '#FFFF33' pt 5 ps 1 lt 1 lw 2\n";
    g << "set style line 7 lc rgb '#A65628' pt 7 ps 1 lt 1 lw 2\n";
    g << "set style line 8 lc rgb '#F781BF' pt 8 ps 1 lt 1 lw 2\n";
    g << "set palette maxcolors 8\n";
    g << "set palette defined ( 0 '#E41A1C', 1 '#377EB8', 2 '#4DAF4A', 3 '#984EA3', 4 '#FF7F00', 5 '#FFFF33', 6 '#A65628', 7 '#F781BF')\n";
    g << "set style line 11 lc rgb '#808080' lt 1 lw 3\n";
    g << "set border 0 back ls 11\n";
    g << "set tics out nomirror\n";
    g << "set style line 12 lc rgb '#808080' lt 0 lw 1\n";
    g << "set grid back ls 12\n";
    g << "set terminal pdfcairo enhanced color dashed font 'Alegreya, 14' rounded size 16 cm, 9.6 cm\n";
  }
  else {
    g << "cd '..'\n";
    g << "load 'stylefile.gp'\n";
  }
}

void Plotter::test(){
  std::cout << "Test Plotter" << std::endl;
  Gnuplot gp;
  setPlotStyle(gp, 0);
  gp << "set output 'test.pdf'\n";
  gp << "set xlabel 'x axis label'\n";
  gp << "set ylabel 'y axis label'\n";
  gp << "set key top right\n";
  gp << "plot [-pi/2:pi] cos(x),-(sin(x) > sin(x+1) ? sin(x) : sin(x+1)))\n";
  gp << "set terminal pop\n";
  gp << "set output\n";
  gp << "replot\n";
}

void Plotter::plot(){
  ;
}

void Plotter::plotPsi(){
  doublePairVec x_psi;
  for(int j = 0; j < system.wavefunctions[0].grid.nPoint; ++j){
    x_psi.push_back(std::make_pair(system.wavefunctions[0].grid.x(j), system.wavefunctions[0].getAbs()(j)));
  }
  Gnuplot gp;
  setPlotStyle(gp, 0);
  gp << "set output 'psi.pdf'\n";
  gp << "set xlabel 'x'\n";
  gp << "set ylabel 'psi'\n";
  gp << "set key top right\n";
  gp << "plot '-' with lines\n";
  gp.send(x_psi);
  gp << "set terminal pop\n";
  gp << "set output\n";
  gp << "replot\n";
}

void Plotter::animate(int nSteps, double stepSize, int evolveOrder, int updateRate){
  Gnuplot gp;
  gp << "set xlabel 'x'\n";
  gp << "set key top right\n";
  gp << "set yrange [-0.5:0.5]\n";
//  gp << "set log y\n";
  for (int j = 0; j < nSteps; ++j) {
    system.evolveAll(stepSize, evolveOrder);
    if (j % updateRate == 0){
      system.log(j*stepSize);
      std::vector<double> xGrid(system.wavefunctions[0].grid.nPoint);
      std::vector<double> psiAbs(system.wavefunctions[0].grid.nPoint);
      std::vector<double> psiReal(system.wavefunctions[0].grid.nPoint);
      std::vector<double> psiImag(system.wavefunctions[0].grid.nPoint);
      VectorXd::Map(&xGrid[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].grid.x;
      VectorXd::Map(&psiAbs[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getAbs();
      VectorXd::Map(&psiReal[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getReal();
      VectorXd::Map(&psiImag[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getImag();
      gp << "plot '-' with lines title 'Norm', '-' with lines title 'Real', '-' with lines title 'Imaginary'\n";
      gp.send(boost::make_tuple(xGrid, psiAbs));
      gp.send(boost::make_tuple(xGrid, psiReal));
      gp.send(boost::make_tuple(xGrid, psiImag));
      gp.flush();
//    pause(pauseTime);
    }
  }
}

void Plotter::animatePsi(int nSteps, double stepSize, int evolveOrder) {
  Gnuplot gp;
  gp << "set xlabel 'x'\n";
//  gp << "set ylabel 'psi'\n";
  gp << "set key top right\n";
  for (int j = 0; j < nSteps; ++j) {
    system.evolveAll(stepSize, evolveOrder);
    doublePairVec x_psi;
    for(int k = 0; k < system.wavefunctions[0].grid.nPoint; ++k){
      x_psi.push_back(std::make_pair((system.wavefunctions[0].grid.x)(k), (system.wavefunctions[0].getAbs())(k)));
    }
    gp << "plot '-' with lines title 'Abs(Psi)'\n";
    gp.send(x_psi);
    gp.flush();
//    pause(pauseTime);
  }
}