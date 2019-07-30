#include "Wavefunction.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <numeric>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "Grid.h"

typedef std::vector<double> double_vec;
typedef std::vector< std::complex<double> > complex_vec;
typedef std::complex<double> complex;

double integrate_vector(double_vec vect, double a, double b, int n){
  double h = (b-a)/n;
  return (h/48.0)*(17*vect[0]+ 59*vect[1]+43*vect[2]+49*vect[3]+48*std::accumulate(vect.begin()+4,vect.end()-3,0)
      +49*vect[n-3]+43*vect[n-2]+59*vect[n-1]+17*vect[n]);
}

Wavefunction::Wavefunction(const Grid& object, double ReducedMass) : grid(1,0.0,1.0,1.0) {
  grid = object;
  reduced_mass = ReducedMass*amu;
  std::cout << "Number of points on grid: " << grid.n_point << std::endl;
  std::cout << "Reduced mass: " << reduced_mass << std::endl;
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

/// Incomplete

double Wavefunction::Overlap(const Wavefunction& object) {
  double_vec integrand;
  for(int i=0; i<grid.n_point; ++i) {
    integrand.push_back(std::abs(psi[i] * std::conj(object.psi[i])));
  }
  return integrate_vector(integrand, grid.x_min, grid.x_max, grid.n_point);
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