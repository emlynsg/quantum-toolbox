#include "Plotter.h"
#include "System.h"

Plotter::Plotter(System sys)
                 : system(Wavefunction(Grid(1, 0.0, 1.0, 1.0), 1.0), Potential(Grid(1, 0.0, 1.0, 1.0))) {
  system = sys;
}

Plotter::~Plotter(){
  std::cout << "Plotter deleted" << std::endl;
}

void Plotter::test(){
  std::cout << "Test Plotter" << std::endl;
  Eigen::ArrayXf x(16);
  x <<   0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
      0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5;

  Eigen::ArrayXf y(16);
  y <<  65,  79,  80,  68,  77,  81, 100, 102,
      105, 111, 120, 126, 120, 104,  85,  92;

  Madplotlib plt;
  plt.plot(x, y);
  plt.show();
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

void Plotter::animate(int nSteps, double stepSize, int evolveOrder, int updateRate, bool logY = false){
  Gnuplot gp;
  gp << "set xlabel 'x'\n";
  gp << "set key top right\n";
  if (logY){
      gp << "set yrange [*:0.5]\n";
      gp << "set log y\n";
      gp << "set format y '10^{%L}'\n";
      for (int j = 0; j < nSteps; ++j) {
        system.evolveAllStep(stepSize, evolveOrder);
        if (j % updateRate == 0){
          system.log(j*stepSize);
          std::vector<double> xGrid(system.wavefunctions[0].grid.nPoint);
          std::vector<double> psiAbs(system.wavefunctions[0].grid.nPoint);
          VectorXd::Map(&xGrid[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].grid.x;
          VectorXd::Map(&psiAbs[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getAbs();
          gp << "plot '-' with lines title 'Norm'\n";
          gp.send(boost::make_tuple(xGrid, psiAbs));
          gp.flush();
      }
    }
  }
  else {
    gp << "set yrange [-0.5:0.5]\n";
    for (int j = 0; j < nSteps; ++j) {
      system.evolveAllStep(stepSize, evolveOrder);
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
}

void Plotter::animateCC(int nSteps, int updateRate, bool logY = false, bool k = false, bool cc = false){
  Gnuplot gp;
  gp << "set xlabel 'x'\n";
  gp << "set key top right\n";
  if (k){
    gp << "set yrange [-4.0:4.0]\n";
    gp << "set xrange [-5.0:5.0]\n";
    for (int j = 0; j < nSteps; ++j) {
      system.evolveCCStep();
      if (j % updateRate == 0){
        system.log(j*system.timeStep);
        gp << "set title 't="+tostring(j*system.timeStep)+"'\n";
        std::vector<double> kGrid(system.wavefunctions[0].grid.nPoint);
        std::vector<double> psiKAbs(system.wavefunctions[0].grid.nPoint);
        std::vector<double> psiKReal(system.wavefunctions[0].grid.nPoint);
        std::vector<double> psiKImag(system.wavefunctions[0].grid.nPoint);
        VectorXd::Map(&kGrid[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].grid.k;
        VectorXd::Map(&psiKAbs[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getKAbs();
        VectorXd::Map(&psiKReal[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getKReal();
        VectorXd::Map(&psiKImag[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getKImag();
        gp << "plot '-' with lines title 'Norm', '-' with lines title 'Real', '-' with lines title 'Imaginary'\n";
        gp.send(boost::make_tuple(kGrid, psiKAbs));
        gp.send(boost::make_tuple(kGrid, psiKReal));
        gp.send(boost::make_tuple(kGrid, psiKImag));
        gp.flush();
      }
    }
    system.updateFromCC();
  }
  else if (logY){
    gp << "set yrange [*:0.5]\n";
    gp << "set log y\n";
    gp << "set format y '10^{%L}'\n";
    for (int j = 0; j < nSteps; ++j) {
      system.evolveCCStep();
      if (j % updateRate == 0){
        system.log(j*system.timeStep);
        gp << "set title 't="+tostring(j*system.timeStep)+"'\n";
        std::vector<double> xGrid(system.wavefunctions[0].grid.nPoint);
        std::vector<double> psiAbs(system.wavefunctions[0].grid.nPoint);
        VectorXd::Map(&xGrid[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].grid.x;
        VectorXd::Map(&psiAbs[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getAbs();
        gp << "plot '-' with lines title 'Norm'\n";
        gp.send(boost::make_tuple(xGrid, psiAbs));
        gp.flush();
      }
    }
    system.updateFromCC();
  }
  else if (cc){
    gp << "set yrange [-0.5:0.5]\n";
//    gp << "set yrange [0.000000000000001:10]\n";
//    gp << "set log y\n";
//    gp << "set format y '10^{%L}'\n";
    for (int j = 0; j < nSteps; ++j) {
      system.evolveCCStep();
      if (j % updateRate == 0){
        system.log(j*system.timeStep);
//        gp << "set title 't="+tostring(j*system.timeStep)+"'\n";
        gp << "set title 'epsilon="+tostring(system.wavefunctions[1].epsilon)+" t="+tostring(j*system.timeStep)+"'\n";
        std::vector<double> xGrid(system.wavefunctions[0].grid.nPoint);
        std::vector<double> pot(system.wavefunctions[0].grid.nPoint);
        std::vector<double> psiAbs(system.wavefunctions[0].grid.nPoint);
        std::vector<double> psiReal(system.wavefunctions[0].grid.nPoint);
        std::vector<double> psiImag(system.wavefunctions[0].grid.nPoint);
        std::vector<double> psiAbs2(system.wavefunctions[0].grid.nPoint);
        std::vector<double> psiReal2(system.wavefunctions[0].grid.nPoint);
        std::vector<double> psiImag2(system.wavefunctions[0].grid.nPoint);
        VectorXd::Map(&xGrid[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].grid.x;
        VectorXd::Map(&psiAbs[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getAbs();
        VectorXd::Map(&psiReal[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getReal();
        VectorXd::Map(&psiImag[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getImag();
        VectorXd::Map(&psiAbs2[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[1].getAbs();
        VectorXd::Map(&psiReal2[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[1].getReal();
        VectorXd::Map(&psiImag2[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[1].getImag();
        VectorXd::Map(&pot[0], system.wavefunctions[0].grid.nPoint) = system.potentials[0].V.real().matrix().normalized().array();
//        gp << "plot '-' with lines title 'Norm', '-' with lines title 'Real', '-' with lines title 'Imaginary', '-' with lines title 'Norm', '-' with lines title 'Real', '-' with lines title 'Imaginary'\n";
        gp << "plot '-' with lines title 'Ground', '-' with lines title 'Excited', '-' with lines title 'Potential'\n";
//        gp << "plot '-' with lines title 'G', '-' with lines title 'G_{Re}', '-' with lines title 'G_{Im}', '-' with lines title 'Ex', '-' with lines title 'Ex_{Re}', '-' with lines title 'Ex_{Im}', '-' with lines title 'Potential'\n";
        gp.send(boost::make_tuple(xGrid, psiAbs));
//        gp.send(boost::make_tuple(xGrid, psiReal));
//        gp.send(boost::make_tuple(xGrid, psiImag));
        gp.send(boost::make_tuple(xGrid, psiAbs2));
//        gp.send(boost::make_tuple(xGrid, psiReal2));
//        gp.send(boost::make_tuple(xGrid, psiImag2));
        gp.send(boost::make_tuple(xGrid, pot));
        gp.flush();
      }
    }
    system.updateFromCC();
  }
  else {
    gp << "set yrange [-0.5:0.5]\n";
    for (int j = 0; j < nSteps; ++j) {
      system.evolveCCStep();
      if (j % updateRate == 0){
        system.log(j*system.timeStep);
        gp << "set title 't="+tostring(j*system.timeStep)+"'\n";
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
      }
    }
    system.updateFromCC();
  }
}

void Plotter::animatePsi(int nSteps, double stepSize, int evolveOrder) {
  Gnuplot gp;
  gp << "set xlabel 'x'\n";
//  gp << "set ylabel 'psi'\n";
  gp << "set key top right\n";
  for (int j = 0; j < nSteps; ++j) {
    system.evolveAllStep(stepSize, evolveOrder);
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