#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <vector>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_complex.h>

#include "Grid.h"

class Potential {
 public:
  /// Objects ///
  Grid grid;
  /// Potential ///
  complexVec V;
  /// Functions ///
  explicit Potential(const Grid &object);
  ~Potential();
  void test();
  doubleVec getReal();
  doubleVec getImag();
  doubleVec getAbs();
  void initZero();
  void initConstantInRegion(const double &c, const double &xmin, const double &xmax);
  void addConstant(const double &c, const double &xmin, const double &xmax);
  void addParabolic(const double &xCenter, const double &c);
  void addQuartic(const double &xCenter, const double &c);
  void addGaussian(const double &xCenter, const double &height, const double &sigma);
  void addWoodsSaxon(const double &xCenter, const double &height, const double &xSize, const double &diffuseness);
  void addCoulomb(const double &Z1Z2, const double &xCenter, const double &xSize);
};

#endif //POTENTIAL_H
