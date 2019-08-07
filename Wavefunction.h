
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

typedef std::complex<double> complex;
typedef std::vector< complex > complex_vec;

class Wavefunction {

 public:

  /// Constants ///
  double hbarc=197.3; /// MeV fm
  double amu=931.5;   /// MeV/c^2
  complex i = complex (0, 1);
  /// Parameters ///
  double reduced_mass;
  /// Objects ///
  Grid grid;
  /// Psi ///
  complex_vec psi;
  complex_vec psi_k;

  /// Functions ///
  Wavefunction(const Grid& object, const double& ReducedMass);
  ~Wavefunction();
  void TestFcn();
  double Overlap(const Wavefunction& object);
  double Norm();
  void Normalise();
  double NormInRegion(const double& xmin, const double& xmax);
  void ComputePsiK();
  void ComputePsi();
  void Init_Zero();
  void Init_Gaussian(const double& mean, const double& sigma);
  void Init_AsymGaussian(const double& mean, const double& sigma);
  void ZeroEdges();
  void Init_Sine(const double& N);
  void Init_Constant();
  void Boost_WaveNumber(const double& WN);
  void Boost_Energy(const double& energy);
  double_vec Get_Real();
  double_vec Get_Imag();
  double_vec Get_Abs();
  double_vec Get_AbsSq();
  double_vec Get_K_Real();
  double_vec Get_K_Imag();
  double_vec Get_K_Abs();
  double_vec Get_K_AbsSq();
  double Get_AvgX();
  void Copy(const Wavefunction& wf);



};

#endif //WAVEFUNCTION_H
