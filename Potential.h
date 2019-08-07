#ifndef POTENTIAL_H
#define POTENTIAL_H

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

typedef std::complex<double> complex;
typedef std::vector< complex > complex_vec;

class Potential {
 public:

  /// Constants ///
  double esq=1.44;    // Electron charge in MeV fm
  /// Objects ///
  Grid grid;
  /// Potential ///
  complex_vec V;

  /// Functions ///
  explicit Potential(const Grid& object);
  ~Potential();
  void TestFcn();
  double_vec Get_Real();
  double_vec Get_Imag();
  double_vec Get_Abs();
  void Init_Zero();
  void Init_ConstantInRegion(const double& c, const double& xmin, const double& xmax);
  void Add_Constant(const double& c, const double& xmin, const double& xmax);
  void Add_Parabolic(const double& xcenter, const double& c);
  void Add_Quartic(const double& xcenter, const double& c);
  void Add_Gaussian(const double& xcenter, const double& height, const double& sigma);
  void Add_WoodsSaxon(const double& xcenter, const double& height, const double& xsize, const double& diffuseness);
  void Add_CoulombSphere(const double& z1z2, const double& xcenter, const double& xsize);
};

#endif //POTENTIAL_H
