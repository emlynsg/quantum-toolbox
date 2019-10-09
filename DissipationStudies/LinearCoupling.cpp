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
#include <chrono>

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

  // Linear scheme defined by
  // N: Number of channels
  // n \in {1,...,N}: Current channel
  // V_{n->n+1}: n*VC Coupling strength
  // Q_n: n*2 Relative energy

  int N = 2;
  double F = 2.0; // Coupling height
  double baseEpsilon = 2.0; // Standard multiplier for epsilon
  std::vector<double> ns;
  std::vector<double> vcs;
  std::vector<double> epsilons;
  for (int j = 0; j < N; ++j) {
    ns.push_back(j);
    vcs.push_back(j*F);
    epsilons.push_back(2.0*baseEpsilon);
  }
  double mu = 1.0;
  double V0 = 100.0;
  double sigma = 5.0;
  double sigmaF = 5.0;
  double couplingCentre = 0.0;


  int time = 4000;
//  double timestep = 1.0;
  double timestep = 0.1;
  //  double timestep = 0.01;
  unsigned int sizeN = 4095;
//  unsigned int sizeN = 8191;
//  unsigned int sizeN = 16383;
  double xmin = -400.0;
  double xmax = 400.0;
  double kscale = 1.0;

  Grid grid(sizeN, xmin, xmax, kscale);

//#pragma omp parallel for
  Wavefunction ground(grid, mu);
  ground.initGaussian(-80.0, 10.0);
  ground.boostEnergy(V0);

  Potential U1(grid);
  U1.initZero();
  U1.addGaussian(0.0, V0, sigma);

  System system(ground, U1);
  for (int j = 1; j < N; ++j) {
    system.addZeroWavefunction(mu, epsilons[j]);
  }
  for (int j = 1; j < N; ++j) {
    system.addGaussianPotential(0.0, V0, sigma, j, j);
  }
  // Nearest neighbour coupling
  for (int j = 1; j < N; ++j) {
    system.addGaussianPotential(couplingCentre, vcs[j], sigma, j, j-1);
    system.addGaussianPotential(couplingCentre, vcs[j], sigma, j-1, j);
  }

  // CC Evolution
  system.initCC(timestep);
//  system.evolveCC(int(time/timestep));
//  system.updateFromCC();
//
    Plotter plot(system);
    plot.animateCC(int(time/timestep), 100, false, false, true);

  dArray T = system.getTransmission();

  std::ofstream Out("LinearCouplingScheme.csv");
  Out << "E" << "," << "T\n";
  for (int j = 1+ground.grid.nPoint/2; j < ground.grid.nPoint; ++j) {
    Out << system.wavefunctions[0].E(j) << "," << T(j) << "\n";
  }
}