#include "Potential.h"

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


Potential::Potential(const Grid& object) : grid(1,0.0,1.0,1.0) {
  grid = object;
  V.resize(grid.n_point, 0.0);
}

Potential::~Potential() {
  complex_vec().swap(V);
  std::cout << "Potential deleted" << std::endl;
}

void Potential::TestFcn() {
  std::cout << "Test Potential" << std::endl;
}

/// Getters

double_vec Potential::Get_Real(){
  double_vec ret_val(grid.n_point);
  std::transform(V.begin(), V.end(), ret_val.begin(), [](auto& elt){return elt.real();});
  return ret_val;
}

double_vec Potential::Get_Imag(){
  double_vec ret_val(grid.n_point);
  std::transform(V.begin(), V.end(), ret_val.begin(), [](auto& elt){return elt.imag();});
  return ret_val;
}

double_vec Potential::Get_Abs(){
  double_vec ret_val(grid.n_point);
  std::transform(V.begin(), V.end(), ret_val.begin(), [](auto& elt){return std::abs(elt);});
  return ret_val;
}

void Potential::Init_Zero(){
  std::fill(V.begin(), V.end(), 0.0);
}

void Potential::Init_ConstantInRegion(const double& c, const double& xmin, const double& xmax){
  for(int i=0; i<grid.n_point; ++i) {
    if(grid.x[i]>=xmin and grid.x[i] <= xmax) {
      V[i] = c;
    }
  }
}

/// Add to potential

void Potential::Add_Constant(const double& c, const double& xmin, const double& xmax){
  for(int i=0; i<grid.n_point; ++i) {
    if(grid.x[i]>=xmin and grid.x[i] <= xmax) {
      V[i] = V[i]+c;
    }
  }
}

void Potential::Add_Parabolic(const double& xcenter, const double& c){
  std::transform(V.begin(), V.end(), grid.x.begin(), V.begin(), [c, xcenter](auto& v, auto& x){return v+(c*pow(x-xcenter,2.0));});
}

void Potential::Add_Quartic(const double& xcenter, const double& c){
  std::transform(V.begin(), V.end(), grid.x.begin(), V.begin(), [c, xcenter](auto& v, auto& x){return v+(c*pow(x-xcenter,4.0));});
}

void Potential::Add_Gaussian(const double& xcenter, const double& height, const double& sigma){
  std::transform(V.begin(), V.end(), grid.x.begin(), V.begin(), [xcenter, height, sigma](auto& v, auto& x){return v+exp(-pow(x-xcenter,2.0)/(2.0*sigma*sigma));});
}

void Potential::Add_WoodsSaxon(const double& xcenter, const double& height, const double& xsize, const double& diffuseness){
  std::transform(V.begin(), V.end(), grid.x.begin(), V.begin(), [xcenter, height, xsize, diffuseness](auto& v, auto& x){return v+(1+exp((abs(x-xcenter)-xsize)/diffuseness));});
}

void Potential::Add_CoulombSphere(const double& z1z2, const double& xcenter, const double& xsize){
  for(int i=0; i<grid.n_point; ++i) {
    if(grid.x[i]-xcenter < xsize) {
      V[i] = V[i]+z1z2*esq*(3.0-pow(std::abs(grid.x[i]-xcenter)/xsize,2.0)/(2.0*xsize));
    }
    else{
      V[i] = V[i] +z1z2*esq/(std::abs(grid.x[i]-xcenter));
    }
  }
}

