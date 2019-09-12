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

#include "Grid.h"
#include "Wavefunction.h"
#include "Potential.h"
#include "System.h"
#include "Plotter.h"

using namespace Eigen;
using namespace std;

int main() {
  std::vector<double> cc_energies;
  std::vector<double> cc_transmissions;
  std::vector<double> cc_norms;

  for (int j = 0; j < 151; ++j) {
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
    cc_energies.push_back(1.0*j);
    cc_transmissions.push_back(system.wavefunctions[0].getNormInRegion(-1.0, 150.0)+system.wavefunctions[1].getNormInRegion(-1.0, 150.0));
    cc_norms.push_back(system.wavefunctions[0].getNorm()+system.wavefunctions[1].getNorm());
  }

  // Output
  std::ofstream CC_Es("Coup_SQ_E.txt");
  std::ofstream CC_Ts("Coup_SQ_T.txt");
  std::ofstream CC_Ns("Coup_SQ_Norm.txt");

  // Push it out
  for (const auto &e : cc_energies) CC_Es << e << "\n";
  for (const auto &e : cc_transmissions) CC_Ts << e << "\n";
  for (const auto &e : cc_norms) CC_Ns << e << "\n";

  return 0;
}
