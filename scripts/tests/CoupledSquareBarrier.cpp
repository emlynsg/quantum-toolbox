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
  Wavefunction ground(grid, ReducedMass);
  ground.initGaussian(-100.0, 10.0);
  ground.boostEnergy(80.0);
  Wavefunction excited(grid, ReducedMass);
  excited.initZero();
  Potential V0(grid);
  V0.initZero();
  V0.addConstant(75.0, -2.0, -1.0);
  Potential VC(grid);
  VC.initZero();
  VC.addConstant(2.0, -2.0, -1.0);
  System system(ground, V0);
  system.addWavefunction(excited);
  system.addPotential(V0, 1, 1);
  system.addPotential(VC, 0, 1);
  system.addPotential(VC, 1, 0);

  // CC Evolution
  system.initCC(0.1);
  Plotter plot(system);
  plot.animateCC(100000, 100);

  return 0;
}
