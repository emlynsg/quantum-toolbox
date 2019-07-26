#ifndef GRID_H
#define GRID_H

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>

typedef std::vector<double> double_vec;
typedef gsl_vector_complex complex_vec;



class Grid {

 public:
  /// Constants ///
  double hbarc=197.3; /// MeV fm
  double amu=931.5;   /// MeV/c^2
  /// Number of steps on the grid ///
  int n_step;
  int n_point;
  /// Grid positions ///
  double x_min;
  double x_max;
  double x_step;
  double_vec x;
  /// Momenta ///
  double k_scale;
  double k_step;
  double k_min;
  double k_max;
  double_vec k;
  /// Energies ///
  double_vec E;

  /// Functions ///
  Grid(int nstep, double xmin, double xmax, double kscale);
  ~Grid();
  void TestFcn();

};

#endif //GRID_H
