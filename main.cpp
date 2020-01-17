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
  unsigned int sizeN = 4095;
  double xmin = -150.0;
  double xmax = 150.0;
  double kscale = 1.0;
  Grid grid(sizeN, xmin, xmax, kscale);
  double ReducedMass = 1.0;
  Wavefunction wavefunction(grid, ReducedMass);
  wavefunction.initGaussian(-100.0, 10.0);
  wavefunction.boostEnergy(80.0);
  Potential potential(grid);
  potential.initZero();
  potential.addConstant(75.0, -2.0, -1.0);
  System system(wavefunction, potential);

  // CC Evolution
  system.initCC(0.1);
  Plotter plot(system);
  plot.test();
}