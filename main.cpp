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

#include "Grid.h"
#include "Wavefunction.h"
#include "Potential.h"
#include "System.h"
#include "Plotter.h"

//#define EIGEN_RUNTIME_NO_MALLOC

using namespace Eigen;
using namespace std;

int main() {
//  Eigen::internal::set_is_malloc_allowed(false);
/// TODO: Fix System so you can add potentials and wavefunctions freely
/// Currently need to add wavefunctions first

  std::vector<double> energies;
  std::vector<double> transmissions;
  std::vector<double> norms;

  for (int j = 50; j < 51; ++j) {
    unsigned int sizeN = 4095;
    double xmin = -150.0;
    double xmax = 150.0;
    double kscale = 1.0;
    Grid grid(sizeN, xmin, xmax, kscale);
    double ReducedMass = 1.0;
    Wavefunction wavefunction(grid, ReducedMass);
    wavefunction.initGaussian(-100.0, 10.0);
    wavefunction.boostEnergy(1.0*j);
    Potential potential(grid);
    potential.initZero();
    potential.addConstant(75.0, -2.0, -1.0);
    System system(wavefunction, potential);
    // CC Evolution
    system.initCC(0.1);
    system.evolveCC(4200);
    system.updateFromCC();
    // Add values to vectors
    energies.push_back(1.0*j);
    transmissions.push_back(system.wavefunctions[0].getNormInRegion(-1.0, 150.0));
    norms.push_back(system.wavefunctions[0].getNorm());
    Plotter plot(system);
    plot.animateCC(1, 1);
  }

  // Output
  std::ofstream Es("SQ_E.txt");
  std::ofstream Ts("SQ_T.txt");
  std::ofstream Ns("SQ_Norm.txt");

  // Push it out
  for (const auto &e : energies) Es << e << "\n";
  for (const auto &e : transmissions) Ts << e << "\n";
  for (const auto &e : norms) Ns << e << "\n";

  return 0;
}

