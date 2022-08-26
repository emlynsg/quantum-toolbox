//
// Created by Emlyn Graham on 9/08/19.
// Test for the Quantum Toolbox library
//
#include <vector>


#include "Grid.h"
#include "Wavefunction.h"
#include "Potential.h"
#include "System.h"
#include "Plotter.h"

using namespace Eigen;
using namespace std;

int main() {
  unsigned int sizeN = 4095;
  double xmin = -150.0;
  double xmax = 150.0;
  double kscale = 1.0;
  Grid grid(sizeN, xmin, xmax, kscale);
  double ReducedMass = 1.0;
  Wavefunction wavefunction(grid, ReducedMass);
  wavefunction.initConstantInRegion(-150.0, -2.0);
  wavefunction.boostEnergy(80.0);
  Potential potential(grid);
  potential.initZero();
  potential.addConstant(75.0, -2.0, -1.0);
  System system(wavefunction, potential);

  // Standard Evolution
  Plotter plot(system);
  plot.animate(100000, 0.1, 20, 100);

  // CC Evolution
  system.initCC(0.1);
  Plotter plot(system);
  plot.animateCC(100000, 100);

  return 0;
}
