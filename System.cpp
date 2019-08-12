#include "System.h"

System::System(Wavefunction &wf, Potential &pot) {
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
  double A = pow(HBARC/wavefunctions[index].grid.xStep,2.0)/(2.0*wavefunctions[index].reducedMass);
  complex B = -1.0*i*timeStep/HBARC;
  complexVec psiPart = wavefunctions[index].psi;
  complexVec psiTemp = wavefunctions[index].psi;
  complexVec psiRotLeft = wavefunctions[index].psi;
  complexVec psiRotRight = wavefunctions[index].psi;
  std::rotate(psiRotLeft.begin(), psiRotLeft.begin()+1, psiRotLeft.end());
  std::rotate(psiRotRight.begin(), psiRotRight.begin() + psiRotRight.size()-1, psiRotRight.end());
  for (int order = 1; order < maxOrder+1; ++order) {
    psiPart = vectorAdd(vectorScale(vectorSubtract(vectorSubtract(vectorScale(psiTemp, 2.0)
              , psiRotLeft), psiRotRight),A),vectorMultiply(potMatrix[index][index].V,psiTemp));
    psiPart[0] = 0.0;
    psiPart[-1] = 0.0;
    psiTemp = vectorScale(vectorScale(psiPart,B),1.0/order);
    wavefunctions[index].psi = vectorAdd(wavefunctions[index].psi, psiTemp);
  }
  wavefunctions[index].zeroEdges();
}

void System::evolveAll(double timeStep, int maxOrder){
  for (int j = 1; j < wavefunctions.size(); ++j) {
    evolve(j, timeStep, maxOrder);
  }
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
  ;
}