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
  // Choose gridstep DeltaX and timestep DeltaT
  std::vector<double> DeltaXs = {1.0, 0.1, 0.01};
  std::vector<double> DeltaTs = {1.0, 0.1, 0.01};

  for (auto const& DeltaX : DeltaXs) {
    for (auto const& DeltaT : DeltaTs) {
      unsigned int sizeN = int(0.5+100.0/DeltaX);
      double xmin = -50.0;
      double xmax = 50.0;
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
      system.evolveCC(200);
    }
  }
  return 0;

}

