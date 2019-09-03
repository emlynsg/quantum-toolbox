//
// Created by Emlyn Graham on 9/08/19.
// System class for performing system evolution.
//

//#include "eigen/unsupported/Eigen/CXX11/Tensor"
#include "eigen/unsupported/Eigen/CXX11/Tensor"

#ifndef SYSTEM_H
#define SYSTEM_H

#include "System.h"
#include "Potential.h"
#include "Wavefunction.h"

typedef std::vector<Wavefunction> wavefunctionVec;
typedef std::vector<Potential> potentialVec;
typedef std::vector<std::vector<Potential>> potentialMatrix;

class System {
 public:
  /// Wavefunction properties ///
  wavefunctionVec wavefunctions; // Vector of wavefunctions
  doubleVecVec times; // Times for each
  doubleVecVec energies; // Energies for each
  doubleVecVec norms; // Norm for each
  doubleVecVec averages; // Average X values for each
  /// Potential properties ///
  potentialVec potentials; // Vector of potentials
  intVec potLeft;
  intVec potRight;
  potentialMatrix potMatrix;
  /// Coupled Channels Objects ///
  cdMatrixTensor potentialTensor;
  cdMatrixTensor U; // Diagonalisation of potential operator
  cdVectorTensor expD;
  cdMatrixTensor Udagger;
  cdVectorTensor expP; // Fourier transformed kinetic operator
  cdVectorTensor psiTensor;

  /// Functions ///
  System(Wavefunction wf, Potential pot);
  ~System();
  void test();
  void addWavefunction(Wavefunction &wf);
  void addPotential(Potential &pot, const int &j, const int &k);
  void evolve(int index, double timeStep, int maxOrder);
  void evolveAll(double timeStep, int maxOrder);
  void initCC();
  void evolveCC(double timeStep);
  void log(double time);
  double energy(int index);
  double hamiltonianElement(int indexI, int indexJ);
};

#endif //SYSTEM_H
