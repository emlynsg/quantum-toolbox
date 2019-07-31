
#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <vector>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <numeric>
#include <functional>
#include <algorithm>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_complex.h>

#include "Grid.h"

typedef std::vector<double> double_vec;
typedef std::complex<double> complex;
typedef std::vector< complex > complex_vec;


class Wavefunction {

 public:

  /// Constants ///
  double hbarc=197.3; /// MeV fm
  double amu=931.5;   /// MeV/c^2
  /// Parameters ///
  double reduced_mass;
  /// Objects ///
  Grid grid;
  /// Psi ///
  complex_vec psi;
  complex_vec psi_k;

  /// Functions ///
  Wavefunction(const Grid& object, double ReducedMass);
  ~Wavefunction();
  void TestFcn();
  double Overlap(Wavefunction& object);
  double Norm();
  void Normalise();
  double NormInRegion(double xmin, double xmax);
  void ComputePsiK();

};

#endif //WAVEFUNCTION_H
