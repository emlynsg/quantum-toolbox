//
// Created by Emlyn Graham on 22/09/19.
// Tests error differences for various time and position steps
//
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <fstream>
#include <cmath>
#include "omp.h"
#include <string>

#include "Grid.h"
#include "Wavefunction.h"
#include "Potential.h"
#include "System.h"
#include "Plotter.h"

//#define EIGEN_RUNTIME_NO_MALLOC
using namespace Eigen;
using namespace std;

// Analytic transmission for 5fm barrier width, 1AMU reduced mass.
// Expression computed using Mathematica
double AnalyticTransmission(double RelativeE, double V0, double width, double mass){
  if(RelativeE == 1.0){
    return HBARC*HBARC/(0.5*V0*width*width*mass+HBARC*HBARC);
  }
  else {
    if(RelativeE < 1.0){
      return 1.0/(1.0+(V0*pow(sinh(sqrt(2.0)*width*sqrt(((V0-RelativeE*V0)*mass)/(HBARC*HBARC))),2.0))/(4.0*RelativeE*(V0-RelativeE*V0)));
    }
    else {
      return 1.0/(1.0+(V0*pow(sin(sqrt(2.0)*width*sqrt(((RelativeE*V0-V0)*mass)/(HBARC*HBARC))),2.0))/(4.0*RelativeE*(RelativeE*V0-V0)));
    }
  }
}

int main() {
  omp_set_num_threads(4);
  Eigen::setNbThreads(4);
//  Eigen::internal::set_is_malloc_allowed(false);
/// TODO: Fix System so you can add potentials and wavefunctions freely
/// Currently need to add wavefunctions first
  int time = 1000;
  double timestep = 0.1;
//  double timestep = 0.01;

//  std::vector<double> Fs = {0.0, 2.0};

  unsigned int sizeN = 2047;
//  unsigned int sizeN = 16383;
  double xmin = -200.0;
  double xmax = 200.0;
  double kscale = 1.0;
  Grid grid(sizeN, xmin, xmax, kscale);
  std::vector<cd> initPsiK;
  double mu = 1.0;
  double V1 = 10.0;
  double V2 = 10.0;
  double sigma1 = 6.0;
  double sigma2 = 6.0;
  double F = 2.0; // coupling potential amplitude
  std::ofstream Out("CoupledWF.csv");

  Wavefunction ground(grid, mu);

  ground.initGaussian(-50.0, 10.0);
  ground.boostEnergy(V1);
  Wavefunction excited(grid, mu);
  excited.initZero();

  Potential U1(grid);
  U1.initZero();
  U1.addGaussian(0.0, V1, 6.0);

  Potential U2(grid);
  U2.initZero();
  U2.addGaussian(0.0, V2, 6.0);

  Potential VC(grid);
  VC.initZero();
  VC.addGaussian(0.0, V2, 6.0);

  System system(ground, U1);
  system.addWavefunction(excited);
  system.addPotential(VC, 0, 1);
  system.addPotential(VC, 1, 0);
  system.addPotential(U2, 1, 1);

  system.updateK();
  dArray groundPsiK_init = system.wavefunctions[0].psiK.abs();
  dArray excitedPsiK_init = system.wavefunctions[1].psiK.abs();

  // CC Evolution
  system.initCC(timestep);
  system.evolveCC(int(time/timestep));
  system.updateFromCC();
//
//  Plotter plot(system);
//  plot.animateCC(int(time/timestep), 100, false, false, true);
  Out << "x" << "," << "InitGround" << "," << "InitExcited" << "," << "FinalGround" << "," << "FinalExcited\n";
  for (int j = 0; j < grid.nPoint; ++j) {
    Out << grid.x(j) << "," << groundPsiK_init(j) << "," << excitedPsiK_init(j) << "," << system.wavefunctions[0].psi.abs()(j) << "," << system.wavefunctions[1].psi.abs()(j) << "\n";
  }


}