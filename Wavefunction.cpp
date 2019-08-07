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

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

struct exponentiate{complex operator()(complex d)const{return std::exp(d);}};

double Sine(const double& x, const double& N, const double& xmin, const double& xmax){
  return sin(N*M_PI*(x-xmin)/(xmax-xmin));
}

double Gaussian(const double& x, const double& X0, const double& Sigma){
  return exp(-pow(x-X0,2.0)/(2*Sigma*Sigma));
}

double AsymGaussian(const double& x, const double& X0, const double& Sigma) {
  return 1.0*(exp(-pow(x-X0-Sigma,2.0)/(2*Sigma*Sigma))-exp(-pow(x-X0+Sigma,2.0)/(2*Sigma*Sigma)));
}

double_vec complex_to_double(const complex_vec& a){
  double_vec b(a.size());
  std::transform(a.begin(), a.end(), b.begin(), [](complex elt){return elt.real();});
  return b;
}

complex_vec double_to_complex(const double_vec& a){
  complex_vec b(a.size());
  std::transform(a.begin(), a.end(), std::back_inserter(b),
      [](double r) { return std::complex<double>(r, 0.0); });
  return b;
}

complex_vec multiply_vecs(const complex_vec& a, const complex_vec& b){
  assert(("Vector lengths don't match", a.size() == b.size()));
  complex_vec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::multiplies<>());
  return c;
}

double_vec multiply_vecs(const double_vec& a, const double_vec& b){
  assert(("Vector lengths don't match", a.size() == b.size()));
  double_vec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::multiplies<>());
  return c;
}

complex_vec exp_vec(const complex_vec& a){
  complex_vec b(a.size());
  std::transform(a.begin(), a.end(), b.begin(), exponentiate());
  return b;
}

complex_vec scale_vec(const complex_vec& a, const complex& b){
  complex_vec c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto& elt){return elt*b;});
  return c;
}

complex_vec scale_vec(const complex_vec& a, const double& b){
  complex_vec c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto& elt){return elt*b;});
  return c;
}

double_vec scale_vec(const double_vec& a, const double& b){
  double_vec c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto& elt){return elt*b;});
  return c;
}

complex_vec scale_vec(const double_vec& a, const complex& b){
  complex_vec c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto& elt){return elt*b;});
  return c;
}


complex_vec add_vecs(const complex_vec& a, const complex_vec& b){
  complex_vec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<>());
  return c;
}

complex_vec add_vecs(const complex_vec& a, const double_vec& b){
  complex_vec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<>());
  return c;
}

complex_vec add_vecs(const double_vec& a, const complex_vec& b){
  complex_vec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<>());
  return c;
}



double_vec complex_vec_to_double_vec_fft(const complex_vec& cvect){
  double_vec dvect;
  dvect.reserve((2*cvect.size()));
  for(auto i : cvect){
    dvect.push_back(i.real());
    dvect.push_back(i.imag());
  }
  return dvect;
}

complex_vec double_vec_to_complex_vec_fft(const double_vec& dvect){
  complex_vec cvect;
  cvect.reserve((dvect.size())/2);
  int range = 0;
  while (range < dvect.size()){
    complex c = complex(dvect[range], dvect[range+1]);
    cvect.push_back(c);
    range = range + 2;
  }
  return cvect;
}

/// Simpson Rule (from Wikipedia, not sure of reference)
/// might need changing later

double simp_integrate_vector(const double_vec& vect, const double& h, const int& n){
  assert(("Integration requires a minimum of 9 points", n > 9));
  return (h/48.0)*(17.0*vect[0]+ 59.0*vect[1]+43.0*vect[2]+49.0*vect[3]
      +48.0*std::accumulate(vect.begin()+4,vect.begin()+(vect.size()-4),0.0)
      +49.0*vect[n-3]+43.0*vect[n-2]+59.0*vect[n-1]+17.0*vect[n]);
}


Wavefunction::Wavefunction(const Grid& object, const double& ReducedMass) : grid(1,0.0,1.0,1.0) {
  grid = object;
  reduced_mass = ReducedMass*amu;
  psi_k.resize(grid.n_point, 0.0);
  psi.resize(grid.n_point, 0.0);
}

Wavefunction::~Wavefunction() {
  complex_vec().swap(psi);
  complex_vec().swap(psi_k);
  std::cout << "Wavefunction deleted" << std::endl;
}

void Wavefunction::TestFcn() {
  std::cout << "Test Wavefunction" << std::endl;
}

double Wavefunction::Norm() {
  return Overlap((*this));
}

void Wavefunction::Normalise() {
  double a = sqrt(Norm());
  psi = scale_vec(psi, 1.0/a);
}

double Wavefunction::NormInRegion(const double& xmin, const double& xmax) {
  double_vec integrand;
  for(int i=0; i<grid.n_point; ++i) {
    if(grid.x[i]>=xmin and grid.x[i] <= xmax) {
      integrand.push_back(norm(psi[i]));
    }
  }
  assert(("No points in this range", integrand.size() > 0));
  return simp_integrate_vector(integrand, grid.x_step, integrand.size());
}

/// Check ordering on PsiK and Psi outputs
/// https://www.gnu.org/software/gsl/doc/html/fft.html#overview-of-complex-data-ffts
/// Note: Definitely wrong as ifft then fft doesn't give the original values


void Wavefunction::ComputePsiK() {
  // Compute input for FFT=
  complex_vec psi_input = multiply_vecs(psi,exp_vec(scale_vec(grid.x, -1.*grid.k_min*i)));
  // Need input as a double array to Fourier Transform
  double_vec dvec = complex_vec_to_double_vec_fft(psi_input);
  complex_vec().swap(psi_input);

  // FFT
  gsl_fft_complex_wavetable * wavetable;
  gsl_fft_complex_workspace * workspace;
  wavetable = gsl_fft_complex_wavetable_alloc(grid.n_point);
  workspace = gsl_fft_complex_workspace_alloc(grid.n_point);
  // Check FFT worked
  int res = gsl_fft_complex_inverse(dvec.data(), 1, grid.n_point, wavetable, workspace);
  if(res != 0){
    std::cout << "FFT Failed" << std::endl;
  }
  gsl_fft_complex_wavetable_free(wavetable);
  gsl_fft_complex_workspace_free(workspace);

  //std::rotate(dvec.rbegin(), dvec.rbegin()+int(dvec.size()/2), dvec.rend());
  // Convert back
  psi_k = scale_vec(multiply_vecs(double_vec_to_complex_vec_fft(dvec),exp_vec(scale_vec(grid.k,-1.*grid.x_min*i))), grid.x_step/(sqrt(2.*M_PI)));
  double_vec().swap(dvec);
  // Maybe need to shift entries??
  //std::rotate(psi_k.begin(), psi_k.begin()+int(psi_k.size()/2), psi_k.end());
}

void Wavefunction::ComputePsi() {
  // Compute input for FFT
  complex_vec psi_input = scale_vec(multiply_vecs(psi_k,exp_vec(scale_vec(grid.k, grid.k_min*i)))
                          ,(sqrt(2.*M_PI)/grid.x_step));
  // Need input as a double array to Fourier Transform
  double_vec dvec = complex_vec_to_double_vec_fft(psi_input);
  complex_vec().swap(psi_input);

  //std::rotate(dvec.rbegin(), dvec.rbegin()+int(dvec.size()/2), dvec.rend());

  // FFT
  gsl_fft_complex_wavetable * wavetable;
  gsl_fft_complex_workspace * workspace;
  wavetable = gsl_fft_complex_wavetable_alloc(grid.n_point);
  workspace = gsl_fft_complex_workspace_alloc(grid.n_point);
  // Check FFT worked
  int res = gsl_fft_complex_forward(dvec.data(), 1, grid.n_point, wavetable, workspace);
  if(res != 0){
    std::cout << "FFT Failed" << std::endl;
  }
  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);

  //std::rotate(dvec.begin(), dvec.begin()+int(dvec.size()/2), dvec.end());
  // Convert back
  psi = multiply_vecs(double_vec_to_complex_vec_fft(dvec),exp_vec(scale_vec(grid.x,1.*grid.k_min*i)));
  double_vec().swap(dvec);

  //std::rotate(psi.begin(), psi.begin()+int(psi.size()/2), psi.end());
}

void Wavefunction::Init_Zero(){
  std::fill(psi.begin(), psi.end(), 0.0);
}

void Wavefunction::Init_Gaussian(const double& mean, const double& sigma){
  std::transform(grid.x.begin(), grid.x.end(), psi.begin(), [mean, sigma](auto& elt){return Gaussian(elt, mean, sigma);});
  ZeroEdges();
  Normalise();
}

void Wavefunction::Init_AsymGaussian(const double& mean, const double& sigma){
  std::transform(grid.x.begin(), grid.x.end(), psi.begin(), [mean, sigma](auto& elt){return AsymGaussian(elt, mean, sigma);});
  ZeroEdges();
  Normalise();
}

void Wavefunction::ZeroEdges(){
  psi[0] = complex(0.0,0.0);
  psi.back() = complex(0.0,0.0);
}

void Wavefunction::Init_Sine(const double& N){
  std::transform(grid.x.begin(), grid.x.end(), psi.begin(), [N, this](auto& elt){return Sine(elt, N, grid.x_min, grid.x_max);});
  Normalise();
}

void Wavefunction::Init_Constant(){
  std::fill(psi.begin(), psi.end(), 1.0);
  ZeroEdges();
  Normalise();
}

void Wavefunction::Boost_WaveNumber(const double& WN){
  std::transform(grid.x.begin(), grid.x.end(), psi.begin(), psi.begin(), [WN, this](auto& x, auto& p){return p*exp(i*WN*x);});
  Normalise();
}

void Wavefunction::Boost_Energy(const double& energy){
  double WN = sgn(energy)*sqrt(2*reduced_mass*std::abs(energy))/hbarc;
  Boost_WaveNumber(WN);
}

/// Getters

double_vec Wavefunction::Get_Real(){
  double_vec ret_val(grid.n_point);
  std::transform(psi.begin(), psi.end(), ret_val.begin(), [](auto& elt){return elt.real();});
  return ret_val;
}

double_vec Wavefunction::Get_Imag(){
  double_vec ret_val(grid.n_point);
  std::transform(psi.begin(), psi.end(), ret_val.begin(), [](auto& elt){return elt.imag();});
  return ret_val;
}

double_vec Wavefunction::Get_Abs(){
  double_vec ret_val(grid.n_point);
  std::transform(psi.begin(), psi.end(), ret_val.begin(), [](auto& elt){return std::abs(elt);});
  return ret_val;
}

double_vec Wavefunction::Get_AbsSq(){
  double_vec ret_val(grid.n_point);
  std::transform(psi.begin(), psi.end(), ret_val.begin(), [](auto& elt){return pow(std::abs(elt),2.0);});
  return ret_val;
}

double_vec Wavefunction::Get_K_Real(){
  double_vec ret_val(grid.n_point);
  std::transform(psi_k.begin(), psi_k.end(), ret_val.begin(), [](auto& elt){return elt.real();});
  return ret_val;
}

double_vec Wavefunction::Get_K_Imag(){
  double_vec ret_val(grid.n_point);
  std::transform(psi_k.begin(), psi_k.end(), ret_val.begin(), [](auto& elt){return elt.imag();});
  return ret_val;
}

double_vec Wavefunction::Get_K_Abs(){
  double_vec ret_val(grid.n_point);
  std::transform(psi_k.begin(), psi_k.end(), ret_val.begin(), [](auto& elt){return std::abs(elt);});
  return ret_val;
}

double_vec Wavefunction::Get_K_AbsSq(){
  double_vec ret_val(grid.n_point);
  std::transform(psi_k.begin(), psi_k.end(), ret_val.begin(), [](auto& elt){return pow(std::abs(elt),2.0);});
  return ret_val;
}

double Wavefunction::Get_AvgX(){
  double_vec integrand(grid.n_point);
  for(int i=0; i<grid.n_point; ++i) {
    integrand[i] = (std::abs(psi[i]*std::conj(psi[i])));
  }
  double ret_val = simp_integrate_vector(multiply_vecs(integrand,grid.x),grid.x_step,grid.n_point);
  double_vec().swap(integrand);
  return ret_val;
}

void Wavefunction::Copy(const Wavefunction& wf){
  grid = wf.grid;
  reduced_mass = wf.reduced_mass;
  psi.reserve(wf.grid.n_point);
  psi_k.reserve(wf.grid.n_point);
  psi = wf.psi;
  psi_k = wf.psi_k;
}


/// Incomplete

double Wavefunction::Overlap(const Wavefunction& object) {
  double_vec integrand(grid.n_point);
  for(int i=0; i<grid.n_point; ++i) {
    integrand[i] = (std::abs(psi[i]*std::conj(object.psi[i])));
  }
  double ret_val = simp_integrate_vector(integrand, grid.x_step, grid.n_point);
  double_vec().swap(integrand);
  return ret_val;
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