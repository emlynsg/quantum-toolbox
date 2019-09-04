#include "System.h"
#include <eigen/Eigen/Eigenvalues>

System::System(Wavefunction wf, Potential pot) {
  matrixContraction = { Eigen::IndexPair<int>(1, 0) };
  addWavefunction(wf);
  addPotential(pot, 0, 0);
}

System::~System() {
  std::cout << "System deleted" << std::endl;
}

void System::test() {
  std::cout << "Test System" << std::endl;
}

void System::addWavefunction(Wavefunction &wf) {
  wavefunctions.push_back(wf);
  times.push_back(doubleVec());
  energies.push_back(doubleVec());
  norms.push_back(doubleVec());
  averages.push_back(doubleVec());
}

void System::addPotential(Potential &pot, const int &j, const int &k) {
  potentials.push_back(pot);
  potLeft.push_back(j);
  potRight.push_back(k);
  potMatrix.resize(wavefunctions.size());
  for (int j = 0; j < wavefunctions.size(); ++j)
    potMatrix[j].resize(wavefunctions.size(), Potential(wavefunctions[0].grid));
  for (int a = 0; a < wavefunctions.size(); ++a) {
    for (int b = 0; b < wavefunctions.size(); ++b) {
      for (int c = 0; c < potentials.size(); ++c) {
        if (potLeft[c] == a and potRight[c] == b) {
          potMatrix[a][b].copy(potentials[c]);
        }
      }
    }
  }
  potentialTensor = cdMatrixTensor(wavefunctions.size(), wavefunctions.size(), wavefunctions[0].grid.nPoint);
  potentialTensor.setZero();
  for (int l = 0; l < potentials.size(); ++l) {
    cdMatrix potL = potentials[l].V.matrix();
    (potentialTensor.chip(potLeft[l],0)).chip(potRight[l],0) = Matrix_to_Tensor(potL, wavefunctions[0].grid.nPoint);
  }
}

void System::evolve(int index, double timeStep, int maxOrder){
  /// Taylor expansion method
  /// TODO: Figure approach, label usefully.
  double A = pow((HBARC/(1.0*wavefunctions[index].grid.xStep)),2.0)/(2.0*wavefunctions[index].reducedMass);
  cd B = -1.0*i*timeStep/HBARC;
  cdArray psiPart = wavefunctions[index].psi;
  cdArray psiTemp = wavefunctions[index].psi;
  for (int order = 1; order < maxOrder+1; ++order) {
    cdArray psiRotLeft = psiTemp;
    cdArray psiRotRight = psiTemp;
    psiRotLeft(wavefunctions[index].grid.nStep) = psiTemp(0);
    psiRotRight(0) = psiTemp(wavefunctions[index].grid.nStep);
    for (int j = 0; j < wavefunctions[index].grid.nStep; ++j) {
      psiRotLeft(j) = psiTemp(j+1);
      psiRotRight(j+1) = psiTemp(j);
    }
    psiPart = A*(2.0*psiTemp - psiRotLeft - psiRotRight)+ psiTemp*potMatrix[index][index].V;
    psiPart(0) = 0.0;
    psiPart(wavefunctions[index].grid.nStep) = 0.0;
    psiTemp = psiPart*B/order;
    wavefunctions[index].psi += psiTemp;
  }
  wavefunctions[index].zeroEdges();
}

void System::evolveAll(double timeStep, int maxOrder){
  for (int j = 0; j < wavefunctions.size(); ++j) {
    evolve(j, timeStep, maxOrder);
  }
}

void System::initCC(double timeStep) {
  /// Assuming all based on the same grid, with same size
  /// Initialise with the timestep that will be used for evolution
  psiTensor = cdVectorTensor(wavefunctions.size(), wavefunctions[0].grid.nPoint);
  psiTensor.setZero();
  for (int j = 0; j < wavefunctions.size(); ++j){
    cdMatrix psiMatrix = wavefunctions[j].psi.matrix();
    psiTensor.chip(j,0) = Matrix_to_Tensor(psiMatrix, wavefunctions[j].grid.nPoint);
  }
  U = cdMatrixTensor(wavefunctions.size(), wavefunctions.size(), wavefunctions[0].grid.nPoint);
  U.setZero();
  Udagger = cdMatrixTensor(wavefunctions.size(), wavefunctions.size(), wavefunctions[0].grid.nPoint);
  Udagger.setZero();
  cdVectorTensor D = cdVectorTensor(wavefunctions.size(), wavefunctions[0].grid.nPoint);
  D.setZero();
  // Diagonalisation for finding U, Udagger and expD
  for (int j = 0; j < wavefunctions[0].grid.nPoint; ++j){
    cdVectorTensor potChip = potentialTensor.chip(j,2);
    cdMatrix potMat = Tensor_to_Matrix(potChip, wavefunctions.size(), wavefunctions.size());
    ComplexEigenSolver<MatrixXcd> ces;
    ces.compute(potMat);
    U.chip(j,2) = Matrix_to_Tensor(ces.eigenvectors(), unsigned(wavefunctions.size()), unsigned(wavefunctions.size())) ;// U
    D.chip(j,1) = Vector_to_Tensor(ces.eigenvalues(), unsigned(wavefunctions.size())) ;// D
    cdMatrix Uinv = ces.eigenvectors().inverse();
    Udagger.chip(j,2) = Matrix_to_Tensor(Uinv, unsigned(wavefunctions.size()), unsigned(wavefunctions.size())) ;// Udagger
  }
  expD = ((-i*timeStep*0.5)*D).exp();
  expP = cdVectorTensor(wavefunctions.size(), wavefunctions[0].grid.nPoint);
  expP.setZero();
  // Fourier Transform of Laplacian is momentum operator.
  for (int j = 0; j < wavefunctions.size(); ++j){
    cdVector pVec = (exp((-i*timeStep/(2*wavefunctions[0].reducedMass*HBARC*HBARC))*square(wavefunctions[0].grid.k))).matrix();
    expP.chip(j,0) = Vector_to_Tensor(pVec, wavefunctions[0].grid.nPoint);
  }
  potentialOperator = cdMatrixTensor(wavefunctions.size(), wavefunctions.size(), wavefunctions[0].grid.nPoint);
  cdMatrixTensor UexpD = cdMatrixTensor(wavefunctions.size(), wavefunctions.size(), wavefunctions[0].grid.nPoint);
  for (int row = 0; row < wavefunctions.size(); ++row) {
    UexpD.chip(row, 0) = (U.chip(row, 0))*expD;
  }
  cout << UexpD.contract(Udagger, matrixContraction) << endl;
//  potentialOperator = UexpD.contract(Udagger, matrixContraction);
}

void System::evolveCC(){
  potentialOperator.contract(psiTensor, matrixContraction);
}

void System::log(double time){
  for (int j = 0; j < wavefunctions.size(); ++j) {
    times[j].push_back(std::abs(time));
    energies[j].push_back(energy(j));
    norms[j].push_back(wavefunctions[j].getNorm());
    averages[j].push_back(wavefunctions[j].getAvgX());
  }
}

double System::energy(int index){
  double A = pow(HBARC/wavefunctions[index].grid.xStep,2.0)/(2.0*wavefunctions[index].reducedMass);
  cdArray psiRotLeft = wavefunctions[index].psi;
  cdArray psiRotRight = wavefunctions[index].psi;
  psiRotLeft(wavefunctions[index].grid.nStep) = wavefunctions[index].psi(0);
  psiRotRight(0) = wavefunctions[index].psi(wavefunctions[index].grid.nStep);
  for (int j = 0; j < wavefunctions[index].grid.nStep; ++j) {
    psiRotLeft(j) = wavefunctions[index].psi(j+1);
    psiRotRight(j+1) = wavefunctions[index].psi(j);
  }
  cdArray psiOverlap = A*(2.0*wavefunctions[index].psi-psiRotLeft-psiRotRight)+wavefunctions[index].psi*potMatrix[index][index].V;
  psiOverlap(0) = 0.0;
  psiOverlap(wavefunctions[index].grid.nStep) = 0.0;
  dArray integrand = abs(wavefunctions[index].psi*(psiOverlap.conjugate()));
  double returnValue = vectorTrapezoidIntegrate(integrand, wavefunctions[index].grid.xStep, wavefunctions[index].grid.nStep);
  return returnValue;
}

double System::hamiltonianElement(int indexI, int indexJ){
  double A = pow(HBARC/wavefunctions[indexI].grid.xStep,2.0)/(2.0*wavefunctions[indexI].reducedMass);
  cdArray psiRotLeft = wavefunctions[indexI].psi;
  cdArray psiRotRight = wavefunctions[indexI].psi;
  psiRotLeft(wavefunctions[indexI].grid.nStep) = wavefunctions[indexI].psi(0);
  psiRotRight(0) = wavefunctions[indexI].psi(wavefunctions[indexI].grid.nStep);
  for (int j = 0; j < wavefunctions[indexI].grid.nStep; ++j) {
    psiRotLeft(j) = wavefunctions[indexI].psi(j+1);
    psiRotRight(j+1) = wavefunctions[indexI].psi(j);
  }
  cdArray psiOverlap = A*(2.0*wavefunctions[indexI].psi-psiRotLeft-psiRotRight)+wavefunctions[indexI].psi*potMatrix[indexI][indexJ].V;
  psiOverlap(0) = 0.0;
  psiOverlap(wavefunctions[indexI].grid.nStep) = 0.0;
  dArray integrand = abs(wavefunctions[indexI].psi*(psiOverlap.conjugate()));
  double returnValue = vectorTrapezoidIntegrate(integrand, wavefunctions[indexI].grid.xStep, wavefunctions[indexI].grid.nStep);
  return returnValue;
}