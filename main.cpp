//
// Created by Emlyn Graham on 9/08/19.
// Main tests all aspects of the Quantum Toolbox library
//

//Standard includes

#include <vector>

// GSL includes

#include "Grid.h"
#include "Wavefunction.h"
#include "Potential.h"
#include "System.h"
#include "Plotter.h"
#include "eigen/Eigen/Dense"
#include <chrono>

using namespace Eigen;
using namespace std;

std::array<double, 1024> transformArray(std::array<double, 1024> &array){
  std::array<double, 1024> ret;
  std::transform(array.begin(), array.end(), ret.begin(), [](auto &elt) { return 2.5*elt; });
  return ret;
}

void alterArray(std::array<double, 1024> &array, std::array<double, 1024> &out){
  std::transform(array.begin(), array.end(), out.begin(), [](auto &elt) { return 2.5*elt; });
}

std::vector<double> transformVector(std::vector<double> &vec){
  std::vector<double> ret(1024);
  std::transform(vec.begin(), vec.end(), ret.begin(), [](auto &elt) { return 2.5*elt; });
  return ret;
}

void alterVector(std::vector<double> &vec, std::vector<double> &out){
  std::transform(vec.begin(), vec.end(), out.begin(), [](auto &elt) { return 2.5*elt; });
}

ArrayXd transformEigen(ArrayXd &eig){
  ArrayXd ret = 2.5*eig;
  return ret;
}

void alterEigen(ArrayXd &eig, ArrayXd &out){
  out = 2.5*eig;
}

void test(){
  auto starteig = std::chrono::high_resolution_clock::now();
  for (int k = 0; k < 100000; ++k) {
    ArrayXd testeig = ArrayXd::LinSpaced(1, 0, 1023);
    ArrayXd outeig;
    alterEigen(testeig, outeig);
  }
  auto finisheig = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> teig = finisheig - starteig;
  auto starteig2 = std::chrono::high_resolution_clock::now();
  for (int k = 0; k < 100000; ++k) {
    ArrayXd testeig = ArrayXd::LinSpaced(1, 0, 1023);
    ArrayXd outeig;
    outeig = transformEigen(testeig);
  }
  auto finisheig2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> teig2 = finisheig2 - starteig2;
  std::cout << " Eigen inplace: " << teig.count() << " Eigen out-of-place: " << teig2.count() << std::endl;
}

int main() {
  /// Error example: Seems to go wrong when changing psiPart inside the Taylor expansion loop in evolve from System.cpp
  int sizeN = 1023;
  double xmin = -200.0;
  double xmax = 200.0;
  double kscale = 1.0;
  Grid gridObject(sizeN, xmin, xmax, kscale);
  double ReducedMass = 1;
  Wavefunction waveObject(gridObject, ReducedMass);
  waveObject.initConstant();
  Potential pot(gridObject);
  pot.initZero();
  pot.addParabolic(0.0, 20.0);
  System sys(waveObject, pot);
  Plotter plot(sys);
  plot.animate(10, 0.1, 20);


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

