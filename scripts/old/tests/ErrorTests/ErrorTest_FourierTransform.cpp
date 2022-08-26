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

  // Test for FFT error

  std::vector<double> DeltaXs = {0.001, 0.01, 0.1, 1.0};
#pragma omp parallel for
  for (int j = 0; j < DeltaXs.size(); ++j) {
    double DeltaX = DeltaXs[j];
    unsigned int sizeN = int(lround(1000.0/DeltaX));
    double xmin = -500.0;
    double xmax = 500.0;
    double kscale = 1.0;
    Grid grid(sizeN, xmin, xmax, kscale);
    double ReducedMass = 1.0;
    Wavefunction wavefunction(grid, ReducedMass);
    wavefunction.initGaussian(0.0, 20.0);

    // Save initial wavefunction
    cdArray initPsi = wavefunction.psi;

    wavefunction.computePsiK();
    wavefunction.computePsi();

    cdArray finalPsi = wavefunction.psi;
    std::ofstream Out("FFTError_dX_"+tostring(DeltaX)+".csv");
    Out << "x" << "," << "initPsiReal" << "," << "initPsiImag" << "," << "finalPsiReal" << "," << "finalPsiImag\n";
    for (int k = 0; k < finalPsi.size(); ++k) {
      Out << grid.x(k) << "," << initPsi(k).real() << "," << initPsi(k).imag() << "," << finalPsi(k).real() << "," << finalPsi(k).imag() << "\n";
    }
  }

}