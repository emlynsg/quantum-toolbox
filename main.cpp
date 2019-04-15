#include <iostream>
#include <vector>
#include "Grid.h"
#include "Wavefunction.h"
typedef std::vector<double> double_vec;

int main() {

  /// We need to first make some variables for our Grid class, and then put them into the constructor ///

  int sizeN = 10;
  double xmin = 0.0;
  double xmax = 1.0;
  double kscale = 1.0;

  Grid gridObject(sizeN, xmin, xmax, kscale);

  /// Checking that the Grid class obj ect was instantiated properly ///

  gridObject.TestFcn();
  std::cout << gridObject.x[4] << std::endl;

  /// Checking the Wavefunction class object was instantiated properly ///
  double ReducedMass = 1;
  Wavefunction waveObject(gridObject, ReducedMass);
  waveObject.TestFcn();
  for(int i=0; i<waveObject.grid.n_point; ++i) {
    std::cout << waveObject.grid.x[i] << " ";
  }
  std::cout << std::endl;
  for(int i=0; i<waveObject.grid.n_point; ++i) {
    std::cout << waveObject.psi[i] << " ";
  }
  std::cout << std::endl;

  return 0;

}