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

std::array<double, 1024> transformArray(std::array<double, 1024> &array){
  std::array<double, 1024> ret;
  std::transform(array.begin(), array.end(), ret.begin(), [](auto &elt) { return 2.5*elt; });
}

void alterArray(std::array<double, 1024> &array){
  std::transform(array.begin(), array.end(), array.begin(), [](auto &elt) { return 2.5*elt; });
}

std::vector<double> transformVector(std::vector<double> &vec){
  std::vector<double> ret(1024);
  std::transform(vec.begin(), vec.end(), ret.begin(), [](auto &elt) { return 2.5*elt; });
}

void alterVector(std::vector<double> &vec){
  std::transform(vec.begin(), vec.end(), vec.begin(), [](auto &elt) { return 2.5*elt; });
}

void test(){
  for (int k = 0; k < 100; ++k) {
    std::vector<double> testvec(1024);
    std::array<double, 1024> testarray ;
    for (int j = 0; j < 1024; ++j) {
      testvec[j] = j;
      testarray[j] = j;
    }
    testarray = transformArray(testarray);
    testvec = transformVector(testvec);
    alterArray(testarray);
    alterVector(testvec);
  }

}

int main() {
  /// Consider using C++17 for the parallel operations
  /// https://www.bfilipek.com/2018/11/parallel-alg-perf.html


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

