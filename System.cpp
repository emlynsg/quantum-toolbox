#include "System.h"

System::System(Wavefunction wf, Potential pot) {
  addWavefunction(wf);
  addPotential(pot, 0, 0);
}

System::~System() {
  wavefunctionVec().swap(wavefunctions);
  doubleVecVec().swap(times);
  doubleVecVec().swap(energies);
  doubleVecVec().swap(norms);
  doubleVecVec().swap(averages);
  potentialVec().swap(potentials);
  intVec().swap(potLeft);
  intVec().swap(potRight);
  potentialMatrix().swap(potMatrix);
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
  potLeft.push_back(k);
  potRight.push_back(j);
  potMatrix.resize(wavefunctions.size(), potentialVec(wavefunctions.size(), Potential(Grid(1, 0.0, 1.0, 1.0))));
  for (int a = 0; a < wavefunctions.size(); ++a) {
    for (int b = 0; b < wavefunctions.size(); ++b) {
      for (int c = 0; c < potentials.size(); ++c) {
        if ((potLeft[c] == a and potRight[c] == b) or (potLeft[c] == b and potRight[c] == a)) {
          potMatrix[a][b] = potentials[c];
        }
      }
    }
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

void System::initCC() {
  /// Assuming all based on the same grid, with same size
  psiTensor = cdVectorTensor(wavefunctions.size(), wavefunctions[0].grid.nPoint);
  psiTensor.setZero();
  cout << psiTensor << endl;
  for (int j = 0; j < wavefunctions.size(); ++j){
    cdVector psiMatrix = wavefunctions[j].psi.matrix();
    psiTensor.chip(j,0) = Matrix_to_Tensor(psiMatrix, wavefunctions[j].grid.nPoint);
  }
  U = cdMatrixTensor(wavefunctions.size(), wavefunctions.size(), wavefunctions[0].grid.nPoint);
  U.setZero();
  Udagger = cdMatrixTensor(wavefunctions.size(), wavefunctions.size(), wavefunctions[0].grid.nPoint);
  Udagger.setZero();
  expD = cdVectorTensor(wavefunctions.size(), wavefunctions[0].grid.nPoint);
  expD.setZero();
  expP = cdVectorTensor(wavefunctions.size(), wavefunctions[0].grid.nPoint);
  expP.setZero();
  /// Set up U, Udagger and expD by diagonalization
  for (int j = 0; j < wavefunctions[0].grid.nPoint; ++j) {
//    blah;
//    U.chip(j,2) = ;
//    Udagger.chip(j,2) = ;
//    expD.chip(j,1) = ;
  }
}

void System::evolveCC(double timeStep){
  ;
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
  double returnValue = vectorTrapezoidIntegrate(integrand, wavefunctions[index].grid.xStep, wavefunctions[index].grid.nPoint);
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
  double returnValue = vectorTrapezoidIntegrate(integrand, wavefunctions[indexI].grid.xStep, wavefunctions[indexI].grid.nPoint);
  return returnValue;
}