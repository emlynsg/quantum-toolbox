//
// Created by Emlyn Graham on 9/08/19.
// Includes a class for single-particle wavefunctions.
//

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "Grid.h"

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

class Wavefunction {

 public:
  /// Constants ///
  /// Parameters ///
  double reducedMass;
  /// Objects ///
  Grid grid;
  /// Psi ///
  cdArray psi;
  cdArray psiK;

  /// Functions ///
  Wavefunction(const Grid &object, const double &ReducedMass);
  ~Wavefunction();
  void test();
  double overlap(const Wavefunction &object);
  double getNorm();
  void normalise();
  double getNormInRegion(const double &xmin, const double &xmax);
  void computePsiK();
  void computePsi();
  void initZero();
  void initGaussian(const double &mean, const double &sigma);
  void initAsymmGaussian(const double &mean, const double &sigma);
  void zeroEdges();
  void initSine(const double &N);
  void initConstant();
  void initConstantInRegion(const double &xmin, const double &xmax);
  void boostWaveNumber(const double &WN);
  void boostEnergy(const double &energy);
  dArray getReal();
  dArray getImag();
  dArray getAbs();
  dArray getAbsSq();
  dArray getKReal();
  dArray getKImag();
  dArray getKAbs();
  dArray getKAbsSq();
  double getAvgX();
  void copy(const Wavefunction &wf);
};

#endif //WAVEFUNCTION_H
