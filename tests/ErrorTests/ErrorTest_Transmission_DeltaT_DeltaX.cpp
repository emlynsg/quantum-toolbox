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
    return (4.0*RelativeE*(RelativeE-1.0))/(4.0*RelativeE*(RelativeE-1.0)+pow(sin(width*sqrt(2*V0*mass*(RelativeE-1.0)/(HBARC*HBARC))),2.0));
  }
}

int main() {
  omp_set_num_threads(4);
  Eigen::setNbThreads(4);


  // Testing DeltaT differences
  std::vector<double> DeltaTs = {0.01, 0.1, 1.0, 10.0};
#pragma omp parallel for
  for (int k = 0; k < DeltaTs.size(); ++k) {
    double timestep = DeltaTs[k];
    std::ofstream Out("DeltaT_"+tostring(timestep)+".csv");
    double en = 100.0;
    int time = int(730*sqrt(100.0/en));
    unsigned int sizeN = 16384; // Fix detail in x
    double xmin = -600.0;
    double xmax = 600.0;
    double kscale = 1.0;
    Grid grid(sizeN, xmin, xmax, kscale);
    double ReducedMass = 1.0;
    Wavefunction wavefunction(grid, ReducedMass);
    wavefunction.initGaussian(-30.0, 2.5);
    wavefunction.boostEnergy(en);
    wavefunction.computePsiK();

    // Initial psiK
    cdArray initPsiK = wavefunction.psiK;

    Potential potential(grid);
    potential.initZero();
    potential.addConstant(en, -2.5, 2.5);
    System system(wavefunction, potential);
    // CC Evolution
    system.initCC(timestep);
    system.evolveCC(int(time/timestep));
    system.updateFromCC();
    // Output
    Out << "E/V0" << "," << "AbsPsiKi" << "," << "AbsPsiKf" << "," << "T" << "," << "T_An" << "," << "Error" << "," << "RelError" << "\n";
    dArray T = (abs(system.wavefunctions[0].psiK))/(abs(initPsiK));
    for (int j = wavefunction.grid.nPoint/2 - 1; j < wavefunction.grid.nPoint; ++j) {
      Out << system.wavefunctions[0].grid.E(j)/en << "," << initPsiK.abs()(j) << "," << system.wavefunctions[0].psiK.abs()(j) << "," << T(j) << "," << AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU) << "," << abs(T(j) - AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU)) << "," << abs(T(j) - AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU))/AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU) << "\n";
    }
  }
  // Testing DeltaX differences
  std::vector<unsigned int> sizeNs = {1024, 2048, 4096, 8192, 16384};
#pragma omp parallel for
  for (int k = 0; k < sizeNs.size(); ++k) {
    unsigned int sizeN = sizeNs[k] - 1;
    std::ofstream Out("DeltaX_"+tostring(1200.0/sizeN)+".csv");
    double en = 100.0;
    int time = int(730*sqrt(100.0/en));
    double timestep = 0.1;
    double xmin = -600.0;
    double xmax = 600.0;
    double kscale = 1.0;
    Grid grid(sizeN, xmin, xmax, kscale);
    double ReducedMass = 1.0;
    Wavefunction wavefunction(grid, ReducedMass);
    wavefunction.initGaussian(-30.0, 2.5);
    wavefunction.boostEnergy(en);
    wavefunction.computePsiK();

    // Initial psiK
    cdArray initPsiK = wavefunction.psiK;

    Potential potential(grid);
    potential.initZero();
    potential.addConstant(en, -2.5, 2.5);
    System system(wavefunction, potential);
    // CC Evolution
    system.initCC(timestep);
    system.evolveCC(int(time/timestep));
    system.updateFromCC();
    // Output
    Out << "E/V0" << "," << "AbsPsiKi" << "," << "AbsPsiKf" << "," << "T" << "," << "T_An" << "," << "Error" << "," << "RelError" << "\n";
    dArray T = (abs(system.wavefunctions[0].psiK))/(abs(initPsiK));
    for (int j = wavefunction.grid.nPoint/2 - 1; j < wavefunction.grid.nPoint; ++j) {
      Out << system.wavefunctions[0].grid.E(j)/en << "," << initPsiK.abs()(j) << "," << system.wavefunctions[0].psiK.abs()(j) << "," << T(j) << "," << AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU) << "," << abs(T(j) - AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU)) << "," << abs(T(j) - AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU))/AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU) << "\n";
    }
  }
  // Testing total time differences
  std::vector<double> energies = {40.0, 60.0, 80.0, 100.0, 120.0};
#pragma omp parallel for
  for (int k = 0; k < energies.size(); ++k) {
    double en = energies[k];
    int time = int(730*sqrt(100.0/en));
    std::ofstream Out("TotalT_"+tostring(time)+".csv");
    double timestep = 0.1;
    unsigned int sizeN = 2047;
    double xmin = -600.0;
    double xmax = 600.0;
    double kscale = 1.0;
    Grid grid(sizeN, xmin, xmax, kscale);
    double ReducedMass = 1.0;
    Wavefunction wavefunction(grid, ReducedMass);
    wavefunction.initGaussian(-30.0, 2.5);
    wavefunction.boostEnergy(en);
    wavefunction.computePsiK();

    // Initial psiK
    cdArray initPsiK = wavefunction.psiK;

    Potential potential(grid);
    potential.initZero();
    potential.addConstant(en, -2.5, 2.5);
    System system(wavefunction, potential);
    // CC Evolution
    system.initCC(timestep);
    system.evolveCC(int(time/timestep));
    system.updateFromCC();
    // Output
    Out << "E/V0" << "," << "AbsPsiKi" << "," << "AbsPsiKf" << "," << "T" << "," << "T_An" << "," << "Error" << "," << "RelError" << "\n";
    dArray T = (abs(system.wavefunctions[0].psiK))/(abs(initPsiK));
    for (int j = wavefunction.grid.nPoint/2 - 1; j < wavefunction.grid.nPoint; ++j) {
      Out << system.wavefunctions[0].grid.E(j)/en << "," << initPsiK.abs()(j) << "," << system.wavefunctions[0].psiK.abs()(j) << "," << T(j) << "," << AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU) << "," << abs(T(j) - AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU)) << "," << abs(T(j) - AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU))/AnalyticTransmission(system.wavefunctions[0].grid.E(j)/en, en, 5.0, ReducedMass*AMU) << "\n";
    }
  }
}