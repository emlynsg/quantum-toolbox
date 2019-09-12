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
  int barrier = 40;
  double eMax = 6*barrier;
  double V0 = barrier*1.0;

  int nEnergies = 50*eMax/barrier + 1;

  std::vector<double> energies(nEnergies);
  std::vector<double> transmissions;
  std::vector<double> reflections;
  std::vector<double> norms;

  std::generate(energies.begin(), energies.end(), [n = 0] () mutable { return n++; });
  omp_set_num_threads(6);
  OMP_PROC_BIND="FALSE";
  GOMP_CPU_AFFINITY="1 2 3:2";
  #pragma omp parallel for
  for (auto en: energies) {
    // Time checked for E=80 as 4500 steps, base others on this
    int nSteps = int(4500/sqrt(1.0*en/80.0));
    unsigned int sizeN = 4095;
    double xmin = -200.0;
    double xmax = 200.0;
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
    transmissions.push_back(system.wavefunctions[0].getNormInRegion(-1.0, 200.0));
    reflections.push_back(system.wavefunctions[0].getNormInRegion(-200.0, -2.0));
    norms.push_back(system.wavefunctions[0].getNorm());}

  // Output
  std::ofstream Es("SQ_Barrier_Energies_V0_40.txt");
  std::ofstream Ts("SQ_Barrier_Transmissions_V0_40.txt");
  std::ofstream Rs("SQ_Barrier_Reflections_V0_40.txt");
  std::ofstream Ns("SQ_Barrier_Norms_V0_40.txt");

  // Push it out
  for (const auto &e : energies) Es << e << "\n";
  for (const auto &e : transmissions) Ts << e << "\n";
  for (const auto &e : reflections) Rs << e << "\n";
  for (const auto &e : norms) Ns << e << "\n";

  return 0;
}
