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

  std::vector<int> configs = {1, 2, 3, 4};
#pragma omp parallel for
  for (int k = 0; k < configs.size(); ++k) {
    int N = configs[k];
    double F = 4.0; // Coupling height
    double baseEpsilon = 2.0; // Standard multiplier for epsilon
    std::vector<double> ns;
    std::vector<double> vcs;
    std::vector<double> epsilons;
    for (int j = 0; j < N; ++j) {
      ns.push_back(j);
      vcs.push_back(j*F);
      epsilons.push_back(j*baseEpsilon);
    }
    double mu = 1.0;
    double V0 = 100.0;
    double sigma = 5.0;
    double sigmaF = 5.0;
    double couplingCentre = 0.0;


    int time = 8500;
    double timestep = 1.0;
//    double timestep = 0.1;
    //  double timestep = 0.01;
//  unsigned int sizeN = 4095;
//  unsigned int sizeN = 8191;
//    unsigned int sizeN = 16383;
//    unsigned int sizeN = 32767;
//    unsigned int sizeN = 65535;
    unsigned int sizeN = 131071;
    double xmin = -10000.0;
    double xmax = 10000.0;
    double kscale = 1.0;

    Grid grid(sizeN, xmin, xmax, kscale);

    Wavefunction ground(grid, mu);
    ground.initGaussian(-100.0, 2.5);
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

    // CC EvolutionC_
    system.initCC(timestep);

    auto start = std::chrono::high_resolution_clock::now();
    system.evolveCC(int(time/timestep));
    auto finish = std::chrono::high_resolution_clock::now();

    system.updateFromCC();

//    Plotter plot(system);
//    plot.animateCC(int(time/timestep), 100, false, false, true);

    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time for "+tostring(N)+" couplings: " << elapsed.count() << " s\n";

    dArray T = system.getTransmission();

    std::vector<dArray> Rs;
    for (int l = 0; l < N; ++l) {
      Rs.push_back(system.getReflection(l));
    }

    std::ofstream Out("N_"+tostring(N)+".csv");
    Out << "E,T,";
    for (int l = 0; l < N; ++l) {
      Out << "R"+tostring(l)+",";
    }
    Out << "R\n";
    for (int j = 1+ground.grid.nPoint/2; j < ground.grid.nPoint; ++j) {
      Out << system.wavefunctions[0].E(j) << "," << T(j) << ",";
      double tot = 0.0;
      for (int l = 0; l < N; ++l) {
        double Rl = Rs[l](j);
        Out << Rl << ",";
        tot += Rl;
      }
      Out << tot << "\n";
    }
  }
}