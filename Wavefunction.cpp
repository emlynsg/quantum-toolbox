#include "Wavefunction.h"

//# define NDEBUG
# include <assert.h>

#include <iostream>
#include <vector>
#include <cmath>
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
#include <complex.h>

#include "Grid.h"

typedef std::complex<double> complex;
typedef std::vector< complex > complex_vec;

struct exponentiate{complex operator()(complex d)const{return std::exp(d);}};

double_vec complex_to_double(const complex_vec& a){
  double_vec b;
  std::transform(a.begin(), a.end(), std::back_inserter(b), [](complex elt){return elt.real();});
  return b;
}


complex_vec double_to_complex(const double_vec& a){
  complex_vec b;
  std::transform(a.begin(), a.end(), std::back_inserter(b),
      [](double r) { return std::complex<double>(r, 0.0); });
  return b;
}

complex_vec multiply_vecs(const complex_vec& a, const complex_vec& b){
  assert(("Vector lengths don't match", a.size() == b.size()));
  complex_vec c;
  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(c), std::multiplies<>());
  return c;
}

double_vec multiply_vecs(const double_vec& a, const double_vec& b){
  assert(("Vector lengths don't match", a.size() == b.size()));
  double_vec c;
  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(c), std::multiplies<>());
  return c;
}

complex_vec exp_vec(const complex_vec& a){
  complex_vec b;
  std::transform(a.begin(), a.end(), std::back_inserter(b), exponentiate());
  return b;
}

complex_vec scale_vec(const complex_vec& a, const complex& b){
  complex_vec c;
  std::transform(a.begin(), a.end(), std::back_inserter(c), [b](auto& elt){return elt*b;});
  return c;
}

complex_vec scale_vec(const complex_vec& a, const double& b){
  complex_vec c;
  std::transform(a.begin(), a.end(), std::back_inserter(c), [b](auto& elt){return elt*b;});
  return c;
}

double_vec scale_vec(const double_vec& a, const double& b){
  double_vec c;
  std::transform(a.begin(), a.end(), std::back_inserter(c), [b](auto& elt){return elt*b;});
  return c;
}

complex_vec scale_vec(const double_vec& a, const complex& b){
  complex_vec c;
  std::transform(a.begin(), a.end(), std::back_inserter(c), [b](auto& elt){return elt*b;});
  return c;
}


double_vec complex_vec_to_double_vec_fft(const complex_vec& cvect){
  double_vec dvect;
  for(auto i : cvect){
    dvect.push_back(i.real());
    dvect.push_back(i.imag());
  }
  return dvect;
}

complex_vec double_vec_to_complex_vec_fft(const double_vec& dvect){
  complex_vec cvect;
  int range = 0;
  while (range < dvect.size()){
    complex c = complex (dvect[range], dvect[range+1]);
    cvect.push_back(c);
    range = range + 2;
  }
  return cvect;
}

/// Simpson Rule (from Wikipedia, not sure of reference)
/// might need changing later

double simp_integrate_vector(const double_vec& vect, const double& a, const double& b, const int& n){
  double h = 1.0*(b-a)/(1.0*(n));
  return (h/48.0)*(17.0*vect[0]+ 59.0*vect[1]+43.0*vect[2]+49.0*vect[3]
      +48.0*std::accumulate(vect.begin()+4,vect.begin()+(vect.size()-4),0.0)
      +49.0*vect[n-3]+43.0*vect[n-2]+59.0*vect[n-1]+17.0*vect[n]);
}


Wavefunction::Wavefunction(const Grid& object, const double& ReducedMass) : grid(1,0.0,1.0,1.0) {
  grid = object;
  reduced_mass = ReducedMass*amu;
  for(int i=0; i<grid.n_point; ++i){
    psi.push_back(complex(1.0, 0.0));
    psi_k.push_back(complex(0.0, 0.0));
  }

}

Wavefunction::~Wavefunction() {
  std::cout << "Wavefunction deleted" << std::endl;
}

void Wavefunction::TestFcn() {
  std::cout << "Test Wavefunction" << std::endl;
}

double Wavefunction::Norm() {
  return Overlap(*this);
}

void Wavefunction::Normalise() {
  double a = sqrt(Norm());
  std::transform(psi.begin(), psi.end(), psi.begin(), [a](auto& c){return complex (c.real()/a, c.imag()/a);});
}

double Wavefunction::NormInRegion(const double& xmin, const double& xmax) {
  double_vec integrand;
  for(int i=0; i<grid.n_point; ++i) {
    if(grid.x[i]>xmin and grid.x[i] < xmax){
      integrand.push_back(std::abs(psi[i] * std::conj(psi[i])));
    }
  }
}

/// Check ordering on PsiK output
/// https://www.gnu.org/software/gsl/doc/html/fft.html#overview-of-complex-data-ffts

void Wavefunction::ComputePsiK() {
  /// Compute input for FFT
  complex_vec psi_input = multiply_vecs(psi,exp_vec(scale_vec(grid.x, grid.k_min*i)));
  // Need input as a double array to Fourier Transform
  double_vec dvec = complex_vec_to_double_vec_fft(psi_input);
  int n = int(psi.size());
  double data[2*n];
  std::copy(dvec.begin(), dvec.end(), data);

  // FFT
  gsl_fft_complex_wavetable * wavetable;
  gsl_fft_complex_workspace * workspace;
  wavetable = gsl_fft_complex_wavetable_alloc(n);
  workspace = gsl_fft_complex_workspace_alloc(n);
  // Check FFT worked
  int res = gsl_fft_complex_inverse(data, 1, n, wavetable, workspace);
  if(res != 0){
    std::cout << "FFT Failed" << std::endl;
  }
  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);
}

void Wavefunction::ComputePsi() {

}



/// Incomplete

double Wavefunction::Overlap(Wavefunction& object) {
  double_vec integrand;
  for(int i=0; i<grid.n_point; ++i) {
    integrand.push_back(std::abs(psi[i] * std::conj(object.psi[i])));
  }
  return simp_integrate_vector(integrand, grid.x_min, grid.x_max, grid.n_step);
}
/// All GSL integration seems to need a function as input
/// Strategy: Interpolate points, write this as a function, and then integrate

/*
double Wavefunction::Overlap(const Wavefunction& object) {

  struct spline_parameters{gsl_spline a; gsl_interp_accel b;};

  double grid_array[grid.n_point];
  double integrand[grid.n_point];
  for(int i=0; i<grid.n_point; ++i) {
    grid_array[i] = grid.x[i];
    integrand[i] = std::abs(psi[i] * std::conj(object.psi[i]));
  }
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, grid.n_point);
  gsl_spline_init(spline, grid_array, integrand, grid.n_point);

  double integrand_value = [](double xi, void * parameters) {
    struct spline_parameters *params = (struct spline_parameters *)parameters;
    gsl_spline *splin = &(params->a);
    gsl_interp_accel *ac = &(params->b);
    double integrand_val = gsl_spline_eval(splin, xi, ac);
    return integrand_val;
  };

  gsl_integration_workspace * w
      = gsl_integration_workspace_alloc (1000);

  double result, error;
  gsl_function F;
  struct spline_parameters params;
  params.a = *spline;
  params.b = *acc;
  F.function = &integrand_value;
  F.params = &params;

  gsl_integration_qag(&F, grid.x_min, grid.x_max, 0, 1e-5
                       , grid.n_point, 6, w, &result, &error);

  std::cout << "Result and error: " << result << ", " << error << std::endl;

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  return result;
}


 */