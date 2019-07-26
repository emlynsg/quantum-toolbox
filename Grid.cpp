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
  double_vec xarray(n_point);
  double_vec karray(n_point);
  double_vec earray(n_point);
  for(int i=0; i<n_point; ++i){
    xarray[i] = i*x_step+x_min;
    karray[i] = i*k_step+k_min;
    earray[i] = pow((hbarc*karray[i]),2)/(2.0*amu);
  }
  x = xarray;
  k = karray;
  E = earray;
}

Grid::~Grid() {
  std::cout << "Grid deleted" << std::endl;
}

void Grid::TestFcn() {
  std::cout << "Test Grid" << std::endl;
}