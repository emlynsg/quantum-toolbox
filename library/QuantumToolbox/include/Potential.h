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
  cdArray V;
  /// Functions ///
  explicit Potential(const Grid &object);
  ~Potential();
  void test();
  dArray getReal();
  dArray getImag();
  dArray getAbs();
  void initZero();
  void addConstant(const cd &c, const double &xmin, const double &xmax);
  void addParabolic(const double &xCentre, const cd &c);
  void addQuartic(const double &xCentre, const cd &c);
  void addGaussian(const double &xCentre, const cd &height, const cd &sigma);
  void addWoodsSaxon(const double &xCentre, const double &height, const double &xSize, const double &diffuseness);
  void addCoulomb(const double &Z1Z2, const double &xCentre, const double &xSize);
  void copy(const Potential &pot);
};

#endif //POTENTIAL_H
