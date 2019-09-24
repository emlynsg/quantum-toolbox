//
// Created by Emlyn Graham on 9/08/19.
// Main tests all aspects of the Quantum Toolbox library
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

int main() {
  omp_set_num_threads(4);
  Eigen::setNbThreads(4);
//  Eigen::internal::set_is_malloc_allowed(false);
/// TODO: Fix System so you can add potentials and wavefunctions freely
/// Currently need to add wavefunctions first

  // Test for machine error over time
  double DeltaX = 0.1;
  double DeltaT = 0.1;
  unsigned int sizeN = 2*int(lround(1000.0/DeltaX));
  double xmin = -1000.0;
  double xmax = 1000.0;
  double kscale = 1.0;
  Grid grid(sizeN, xmin, xmax, kscale);
  double ReducedMass = 1.0;
  Wavefunction wavefunction(grid, ReducedMass);
  wavefunction.initGaussian(0.0, 20.0);
  Potential potential(grid);
  potential.initZero();

  // Save initial wavefunction
  dArray initPsi = wavefunction.psi.abs();

  System system(wavefunction, potential);
  // CC Evolution
  system.initCC(DeltaT);
  system.evolveCCStep();
  system.updateFromCC();

  // Save 1-step wavefunction
  dArray Onestep = system.wavefunctions[0].psi.abs();

  system.evolveCC(int(10000/DeltaT));
  system.updateFromCC();

  // Save 1-step wavefunction
  dArray finalPsi = system.wavefunctions[0].psi.abs();


  std::ofstream init("BasicError_InitWF.csv");
  std::ofstream onestep("BasicError_OneStepWF.csv");
  std::ofstream final("BasicError_FinalWF.csv");


  init << "x" << "," << "AbsPsi" << "\n";
  for (int j = 0; j < wavefunction.grid.nPoint; ++j) {
    init << system.wavefunctions[0].grid.x(j) << "," << initPsi(j) << "\n";
  }

  onestep << "x" << "," << "AbsPsi" << "\n";
  for (int j = 0; j < wavefunction.grid.nPoint; ++j) {
    onestep << system.wavefunctions[0].grid.x(j) << "," << Onestep(j) << "\n";
  }

  final << "x" << "," << "AbsPsi" << "\n";
  for (int j = 0; j < wavefunction.grid.nPoint; ++j) {
    final << system.wavefunctions[0].grid.x(j) << "," << finalPsi(j) << "\n";
  }

}