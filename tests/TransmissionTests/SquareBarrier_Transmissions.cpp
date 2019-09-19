//
// Created by Emlyn Graham on 9/08/19.
// Test for the Quantum Toolbox library
//
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <fstream>
#include "omp.h"


#include "Grid.h"
#include "Wavefunction.h"
#include "Potential.h"
#include "System.h"
#include "Plotter.h"

using namespace Eigen;
using namespace std;

int main() {
  omp_set_num_threads(4);
  Eigen::setNbThreads(4);
//  Eigen::internal::set_is_malloc_allowed(false);
/// TODO: Fix System so you can add potentials and wavefunctions freely
/// Currently need to add wavefunctions first

  std::ofstream Out("SQ_Barrier_V0_100.csv");
  std::vector<cd> initPsiK;

  unsigned int sizeN = 16383;
  double xmin = -600.0;
  double xmax = 600.0;
  double kscale = 1.0;
  Grid grid(sizeN, xmin, xmax, kscale);
  double ReducedMass = 1.0;
  Wavefunction wavefunction(grid, ReducedMass);
  wavefunction.initGaussian(-70.0, 2.5);
  wavefunction.boostEnergy(12.0);
  wavefunction.computePsiK();
  for (int k = 0; k < grid.nPoint; ++k) {
    initPsiK.push_back(wavefunction.psiK(k));
  }
  Potential potential(grid);
  potential.initZero();
  potential.addConstant(100.0, -2.5, 2.5);
  System system(wavefunction, potential);
  // CC Evolution
  system.initCC(0.01);
  system.evolveCC(68000);
  system.updateFromCC();

  // Output
  std::vector<double> kGrid(system.wavefunctions[0].grid.nPoint);
  std::vector<double> eGrid(system.wavefunctions[0].grid.nPoint);
  std::vector<double> psiKReal(system.wavefunctions[0].grid.nPoint);
  std::vector<double> psiKImag(system.wavefunctions[0].grid.nPoint);
  VectorXd::Map(&kGrid[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].grid.k;
  VectorXd::Map(&eGrid[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].grid.E;
  VectorXd::Map(&psiKReal[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getKReal();
  VectorXd::Map(&psiKImag[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getKImag();
  Out << "k" << "," << "E" << "," << "RePsiKi" << "," << "ImPsiKi" << "," << "RePsiKf" << "," << "ImPsiKf" << "\n";
  for (int j = 0; j < wavefunction.grid.nPoint; ++j) {
    Out << kGrid[j] << "," << eGrid[j] << "," << initPsiK[j].real() << "," << initPsiK[j].imag() << "," << psiKReal[j] << "," << psiKImag[j] << "\n";
  }

}
