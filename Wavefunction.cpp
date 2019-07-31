#include "Wavefunction.h"

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

#include "Grid.h"

typedef std::vector<double> double_vec;
typedef std::vector< std::complex<double> > complex_vec;
typedef std::complex<double> complex;

#define REAL(z,i) ((z)[2*(i)]) //complex arrays stored as
#define IMAG(z,i) ((z)[2*(i)+1])

/// Simpson Rule (from Wikipedia, not sure of reference)
/// might need changing later

double simp_integrate_vector(double_vec vect, double a, double b, int n){
  double h = 1.0*(b-a)/(1.0*(n));
  return (h/48.0)*(17.0*vect[0]+ 59.0*vect[1]+43.0*vect[2]+49.0*vect[3]
      +48.0*std::accumulate(vect.begin()+4,vect.begin()+(vect.size()-4),0.0)
      +49.0*vect[n-3]+43.0*vect[n-2]+59.0*vect[n-1]+17.0*vect[n]);
}
/*
void complex_vec_to_double_vec(complex_vec cvect, double_vec dvect){
  for(auto i : cvect){
    dvect.push_back((cvect.data()[i]).real());
    dvect.push_back((cvect.data()[i]).imag());
  }
}*/

 /*
gsl_complex_packed_array create_complex_packed_array(complex_vec vect){
  double size[2*vect.size()];
  gsl_complex_packed_array data = size;
  for(int i=0; i<vect.size(); ++i) {
      REAL(data,i) = vect[i].real();
      IMAG(data,i) = vect[i].imag();
  }
  return data;
}

complex_vec read_complex_packed_array(gsl_complex_packed_array array){
  double arr = array;
  int size = (arr.size())/2;
  complex_vec data;
  for(int i=0; i<size; ++i) {
    data.push_back(complex(REAL(array, i), IMAG(array, i)));
  }
  return data;
}
*/

Wavefunction::Wavefunction(const Grid& object, double ReducedMass) : grid(1,0.0,1.0,1.0) {
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

double Wavefunction::NormInRegion(double xmin, double xmax) {
  double_vec integrand;
  for(int i=0; i<grid.n_point; ++i) {
    if(grid.x[i]>xmin and grid.x[i] < xmax){
      integrand.push_back(std::abs(psi[i] * std::conj(psi[i])));
    }
  }
}

void Wavefunction::ComputePsiK() {
  std::cout << psi[0].real() << psi.data()[0] << psi.data()[1] << std::endl;
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