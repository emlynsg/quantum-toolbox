//
// Created by Emlyn Graham on 9/08/19.
// Includes a class for the grid upon which quantum mechanical calculations are performed.
//

#ifndef GRID_H
#define GRID_H

#include "Extras.h"

#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>

#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

class Grid {
 public:
  /// Number of steps on the grid ///
  unsigned int nStep;
  unsigned int nPoint;
  /// Grid positions ///
  double xMin;
  double xMax;
  double xStep;
  dArray x;
  /// Momenta ///
  double kScale;
  double kStep;
  double kMin;
  double kMax;
  dArray k;
  /// Energies ///
  dArray E;

  /// Functions ///
  Grid(unsigned int nstep, double xmin, double xmax, double kscale);
  ~Grid();
  void test();
  void copy(const Grid &gr);
};

#endif //GRID_H
