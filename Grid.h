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

typedef std::vector<double> doubleVec;

class Grid {
 public:
  /// Number of steps on the grid ///
  int nStep;
  int nPoint;
  /// Grid positions ///
  double xMin;
  double xMax;
  double xStep;
  doubleVec x;
  /// Momenta ///
  double kScale;
  double kStep;
  double kMin;
  double kMax;
  doubleVec k;
  /// Energies ///
  doubleVec E;

  /// Functions ///
  Grid(int nstep, double xmin, double xmax, double kscale);
  ~Grid();
  void test();

};

#endif //GRID_H
