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
  std::vector<int> barriers = {10, 40, 75, 100};
  for (auto barrier: barriers) {
    double eMax = 6 * barrier;
    double V0 = barrier * 1.0;
    int nEnergies = int(50 * eMax / barrier + 0.5);
    std::vector<double> cc_energies(nEnergies);
    std::vector<double> cc_transmissions;
    std::vector<double> cc_reflections;
    std::vector<double> cc_norms;
    std::generate(cc_energies.begin(), cc_energies.end(), [&, n = 0]() mutable { return 0.02 * V0 * (n + 1); });
    omp_set_num_threads(3);
#pragma omp parallel for
    for (int j = 0; j < cc_energies.size(); ++j) {
      double en = cc_energies[j];
      // Time checked for E=80 as 4500 steps, base others on this
      int nSteps = int(4500 / sqrt(1.0 * en / 80.0));
      unsigned int sizeN = 4095;
      double xmin = -150.0;
      double xmax = 150.0;
      double kscale = 1.0;
      Grid grid(sizeN, xmin, xmax, kscale);
      double ReducedMass = 1.0;
      Wavefunction wavefunction(grid, ReducedMass);
      wavefunction.initGaussian(-100.0, 10.0);
      wavefunction.boostEnergy(en);
      Potential potential(grid);
      potential.initZero();
      potential.addConstant(V0, -2.0, -1.0);
      System system(wavefunction, potential);
      // CC Evolution
      system.initCC(0.1);
      system.evolveCC(nSteps);
      system.updateFromCC();
      // Add values to vectors
      cc_transmissions.push_back(
          system.wavefunctions[0].getNormInRegion(-1.0, 150.0) + system.wavefunctions[1].getNormInRegion(-1.0, 150.0));
      cc_reflections.push_back(system.wavefunctions[0].getNormInRegion(-150.0, -1.0)
                                   + system.wavefunctions[1].getNormInRegion(-150.0, -1.0));
      cc_norms.push_back(system.wavefunctions[0].getNorm() + system.wavefunctions[1].getNorm());
    }
    // Output
    std::ofstream CC_Es("Coup_SQ_E_V0_" + to_string(barrier) + ".txt");
    std::ofstream CC_Ts("Coup_SQ_T_V0_" + to_string(barrier) + ".txt");
    std::ofstream CC_Rs("Coup_SQ_R_V0_" + to_string(barrier) + ".txt");
    std::ofstream CC_Ns("Coup_SQ_Norm_V0_" + to_string(barrier) + ".txt");
    // Push it out
    for (const auto &e : cc_energies) CC_Es << e << "\n";
    for (const auto &e : cc_transmissions) CC_Ts << e << "\n";
    for (const auto &e : cc_reflections) CC_Rs << e << "\n";
    for (const auto &e : cc_norms) CC_Ns << e << "\n";
  }
  return 0;
}
