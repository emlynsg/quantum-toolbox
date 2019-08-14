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

int main() {

  std::cout << std::endl;

  /// We need to first make some variables for our Grid class, and then put them into the constructor ///

  int sizeN = 1023;
  double xmin = -200.0;
  double xmax = 200.0;
  double kscale = 1.0;

  /// Grid and a pointer to the grid

  Grid gridObject(sizeN, xmin, xmax, kscale);
  Grid *gridPointer;
  gridPointer = &gridObject;

  /// Checking that the Grid class object was instantiated properly ///

  gridObject.test();
  std::cout << "Check grid value: " << gridObject.x[4] << std::endl;


  /// How to print from the pointer

  ///std::cout << (*gridPointer).x[4] << std::endl;

  std::cout << std::endl;
  /// Checking the Wavefunction class object was instantiated properly ///
  double ReducedMass = 1;
  Wavefunction waveObject(gridObject, ReducedMass);
  waveObject.test();
  waveObject.initZero();
  //waveObject.initConstant();
  //waveObject.boostWaveNumber(1.1);
  //waveObject.boostEnergy(1.1);
  //waveObject.initSine(2.5);
  //waveObject.initAsymmGaussian(1.0, 1.0);
  waveObject.initGaussian(0.0, 5.0);
//  doubleVec check = waveObject.getReal();
//  check = waveObject.getImag();
//  check = waveObject.getAbs();
//  check = waveObject.getAbsSq();
//  check = waveObject.getKReal();
//  check = waveObject.getKImag();
//  check = waveObject.getKAbs();
//  check = waveObject.getKAbsSq();
  Wavefunction waveObject2(gridObject, ReducedMass);
  waveObject2.copy(waveObject);
  waveObject.normalise();
//  std::cout << "Norm is: " << waveObject.getNorm() << std::endl;
//  std::cout << "Norm from 0 to 200 is: " << waveObject.getNormInRegion(0.0, 200.0) << std::endl;

  std::cout << "Wavefunction before FFT is: ";
  for (int j = sizeN / 2; j < sizeN / 2 + 5; ++j) {
    std::cout << abs(waveObject.psi[j]) << " ";
  }
  std::cout << std::endl;
  waveObject.computePsiK();
  waveObject.computePsi();
  std::cout << "Wavefunction after FFT and IFFT is: ";
  for (int j = sizeN / 2; j < sizeN / 2 + 5; ++j) {
    std::cout << abs(waveObject.psi[j]) << " ";
  }
  /// Need to fix this Fourier ordering stuff

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

  std::cout << std::endl;

  std::cout << std::endl;
  Potential pot(gridObject);
  pot.test();
  pot.initZero();
//  check = pot.getReal();
//  check = pot.getImag();
//  check = pot.getAbs();
  pot.initConstantInRegion(2.0, -200.0, 200.0);
  std::cout << "Check potential value: " << pot.V[1] << std::endl;
  pot.addConstant(2.0, -200.0, 100.0);
  std::cout << "Check potential value: " << pot.V[1] << std::endl;
  pot.addParabolic(-199.0, 2.0);
  std::cout << "Check potential value: " << pot.V[1] << std::endl;
  pot.addQuartic(-199.0, 2.0);
  std::cout << "Check potential value: " << pot.V[1] << std::endl;
  pot.addGaussian(-198.0, 2.0, 1.0);
  std::cout << "Check potential value: " << pot.V[1] << std::endl;
  pot.addWoodsSaxon(-198.0, 2.0, 1.0, 1.0);
  std::cout << "Check potential value: " << pot.V[1] << std::endl;
  pot.addCoulomb(4.0, -199.5, 0.1);
  std::cout << "Check potential value: " << pot.V[1] << std::endl;

  std::cout << std::endl;

  std::cout << std::endl;

  System sys(waveObject, pot);
  sys.test();
//  sys.evolveAll(10.0, 3);
//  std::cout << "Check Hamiltonian element: " << sys.hamiltonianElement(0,0) << std::endl;



  std::cout << std::endl;

  std::cout << std::endl;

  Plotter plot(sys);
  plot.plotPsi();


  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  return 0;

}