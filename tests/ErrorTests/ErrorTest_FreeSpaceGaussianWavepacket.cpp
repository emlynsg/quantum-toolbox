//
// Created by Emlyn Graham on 9/08/19.
// Test for numerical error changes over propagation scheme
// Test at multiple timesteps and gridsteps
#include <vector>

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

  // Choose gridstep DeltaX and timestep DeltaT
  std::vector<double> DeltaXs = {0.01, 0.1, 1.0};
  std::vector<double> DeltaTs = {0.01, 0.1, 1.0};

#pragma omp parallel for
  for (int k = 0; k < DeltaXs.size(); ++k){
    double DeltaX = DeltaXs[k];
#pragma omp parallel for
    for (int l = 0; l < DeltaTs.size(); ++l) {
      double DeltaT = DeltaTs[l];
      unsigned int sizeN = 2*int(lround(100.0/DeltaX));
      double xmin = -100.0;
      double xmax = 100.0;
      double kscale = 1.0;
      Grid grid(sizeN, xmin, xmax, kscale);
      double ReducedMass = 1.0;
      Wavefunction wavefunction(grid, ReducedMass);
      wavefunction.initGaussian(0.0, 5.0);
      Potential potential(grid);
      potential.initZero();
      System system(wavefunction, potential);
      // CC Evolution
      system.initCC(DeltaT);
      system.evolveCC(int(lround(1000/DeltaT)));
      system.updateFromCC();

      // Output
      std::ofstream Xs("ErrorTest2_x_dt_"+to_string(DeltaT)+"dx_"+to_string(DeltaX)+".txt");
      std::ofstream Norms("ErrorTest2_norm_dt_"+to_string(DeltaT)+"dx_"+to_string(DeltaX)+".txt");
      std::ofstream Reals("ErrorTest2_real_dt_"+to_string(DeltaT)+"dx_"+to_string(DeltaX)+".txt");
      std::ofstream Imags("ErrorTest2_imag_dt_"+to_string(DeltaT)+"dx_"+to_string(DeltaX)+".txt");

      std::vector<double> psiAbs(system.wavefunctions[0].grid.nPoint);
      std::vector<double> psiReal(system.wavefunctions[0].grid.nPoint);
      std::vector<double> psiImag(system.wavefunctions[0].grid.nPoint);
      VectorXd::Map(&psiAbs[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getAbs();
      VectorXd::Map(&psiReal[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getReal();
      VectorXd::Map(&psiImag[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getImag();
      // Push it out
      for (int j = 0; j < system.wavefunctions[0].grid.nPoint; ++j) {
        Xs << system.wavefunctions[0].grid.x(j) << "\n";
        Norms << psiAbs[j] << "\n";
        Reals << psiReal[j] << "\n";
        Imags << psiImag[j] << "\n";
      }
    }
  }

#pragma omp parallel for
  for (int k = 0; k < DeltaXs.size(); ++k){
    double DeltaX = DeltaXs[k];
#pragma omp parallel for
    for (int l = 0; l < DeltaTs.size(); ++l) {
      double DeltaT = DeltaTs[l];
      unsigned int sizeN = 2*int(lround(100.0/DeltaX));
      double xmin = -100.0;
      double xmax = 100.0;
      double kscale = 1.0;
      Grid grid(sizeN, xmin, xmax, kscale);
      double ReducedMass = 1.0;
      Wavefunction wavefunction(grid, ReducedMass);
      wavefunction.initGaussian(0.0, 5.0);
      Potential potential(grid);
      potential.initZero();
      System system(wavefunction, potential);
      // CC Evolution
      system.initCC(DeltaT);
      system.evolveCC(100000);
      system.updateFromCC();

      // Output
      std::ofstream Xs("ErrorTest_x_dt_"+to_string(DeltaT)+"dx_"+to_string(DeltaX)+".txt");
      std::ofstream Norms("ErrorTest_norm_dt_"+to_string(DeltaT)+"dx_"+to_string(DeltaX)+".txt");
      std::ofstream Reals("ErrorTest_real_dt_"+to_string(DeltaT)+"dx_"+to_string(DeltaX)+".txt");
      std::ofstream Imags("ErrorTest_imag_dt_"+to_string(DeltaT)+"dx_"+to_string(DeltaX)+".txt");

      std::vector<double> psiAbs(system.wavefunctions[0].grid.nPoint);
      std::vector<double> psiReal(system.wavefunctions[0].grid.nPoint);
      std::vector<double> psiImag(system.wavefunctions[0].grid.nPoint);
      VectorXd::Map(&psiAbs[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getAbs();
      VectorXd::Map(&psiReal[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getReal();
      VectorXd::Map(&psiImag[0], system.wavefunctions[0].grid.nPoint) = system.wavefunctions[0].getImag();
      // Push it out
      for (int j = 0; j < system.wavefunctions[0].grid.nPoint; ++j) {
        Xs << system.wavefunctions[0].grid.x(j) << "\n";
        Norms << psiAbs[j] << "\n";
        Reals << psiReal[j] << "\n";
        Imags << psiImag[j] << "\n";
      }
    }
  }
}

