#include "Grid.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>

typedef std::vector<double> double_vec;

Grid::Grid(int nstep, double xmin, double xmax, double kscale) {
  n_step = nstep;
  n_point = n_step + 1;
  x_min = xmin;
  x_max = xmax;
  x_step = (x_max - x_min)/n_step;
  k_scale = kscale;
  k_step = 2*M_PI/(n_step*x_step*k_scale);
  k_min = -1*n_step*k_step/2.0;
  k_max = n_step*k_step/2.0;
  x.reserve(n_point);
  k.reserve(n_point);
  E.reserve(n_point);
  for(int i=0; i<n_point; ++i){
    x.push_back(i*x_step+x_min);
    k.push_back(i*k_step+k_min);
    E.push_back(pow((hbarc*k[i]),2.0)/(2.0*amu));
  }
}

Grid::~Grid() {
  double_vec().swap(x);
  double_vec().swap(k);
  double_vec().swap(E);
  std::cout << "Grid deleted" << std::endl;
}

void Grid::TestFcn() {
  std::cout << "Test Grid" << std::endl;
  std::cout << "Number of grid points: " << n_point << std::endl;
}