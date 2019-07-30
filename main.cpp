#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <numeric>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>


#include "Grid.h"
#include "Wavefunction.h"

typedef std::vector<double> double_vec;


int main() {

  /// We need to first make some variables for our Grid class, and then put them into the constructor ///

  int sizeN = 10;
  double xmin = 0.0;
  double xmax = 1.0;
  double kscale = 1.0;

  /// Grid and a pointer to the grid

  Grid gridObject(sizeN, xmin, xmax, kscale);
  Grid *gridPointer;
  gridPointer = &gridObject;

  /// Checking that the Grid class obj ect was instantiated properly ///

  gridObject.TestFcn();
  std::cout << "Check grid value: " << gridObject.x[4] << std::endl;

  /// How to print from the pointer

  ///std::cout << (*gridPointer).x[4] << std::endl;

  std::cout << std::endl;
  /// Checking the Wavefunction class object was instantiated properly ///
  double ReducedMass = 1;
  Wavefunction waveObject(gridObject, ReducedMass);
  waveObject.TestFcn();
  std::cout << "Wavefunction grid is: ";
  for(int i=0; i<waveObject.grid.n_point; ++i) {
    std::cout << waveObject.grid.x[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Wavefunction is: ";
  for(int i=0; i<waveObject.grid.n_point; ++i) {
    std::cout << waveObject.psi[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Overlap is: " << waveObject.Overlap(waveObject) << std::endl;


  /// GSL Matrix Check

  std::cout << "GSL matrix: ";
  int i, j;
  gsl_matrix * m = gsl_matrix_alloc (10, 3);

  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
      gsl_matrix_set (m, i, j, 0.23 + 100*i + j);

  for (i = 0; i < 10; i++)  /* OUT OF RANGE ERROR */
    for (j = 0; j < 3; j++)
      std::cout << "m(" << i << "," << j << ") = " << gsl_matrix_get (m, i, j) << ", ";
  std::cout << std::endl;

  gsl_matrix_free (m);


  return 0;

}