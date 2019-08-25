//
// Created by Emlyn Graham on 9/08/19.
// Main tests all aspects of the Quantum Toolbox library
//

#include <vector>


#include "Grid.h"
#include "Wavefunction.h"
#include "Potential.h"
//#include "System.h"
//#include "Plotter.h"

using namespace Eigen;
using namespace std;

int main() {
  /// Error example: Seems to go wrong when changing psiPart inside the Taylor expansion loop in evolve from System.cpp
  unsigned int sizeN = 1024;
  double xmin = -200.0;
  double xmax = 200.0;
  double kscale = 1.0;
  Grid grid(sizeN, xmin, xmax, kscale);
  double ReducedMass = 1.0;
  Wavefunction wavefunction(grid, ReducedMass);
  wavefunction.initGaussian(0.0, 10.0);
  Potential potential(grid);
  potential.initZero();
  potential.addParabolic(0.0, 20.0);
//  System sys(waveObject, pot);
//  Plotter plot(sys);
//  plot.animate(10, 0.1, 20);


  /*
 /// GSL Matrix Check

 std::cout << "GSL matrix: ";
 int j, k;
 gsl_matrix * m = gsl_matrix_alloc (10, 3);

 for (j = 0; j < 10; j++)
   for (k = 0; k < 3; k++)
     gsl_matrix_set (m, j, k, 0.23 + 100*j + k);

 for (j = 0; j < 10; j++)  // OUT OF RANGE ERROR
   for (k = 0; k < 3; k++)
     std::cout << "m(" << j << "," << k << ") = " << gsl_matrix_get (m, j, k) << ", ";
 std::cout << std::endl;

 gsl_matrix_free (m);

 */

  return 0;

}

