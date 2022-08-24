#include "System.h"
#include <eigen3/unsupported/Eigen/FFT>
#define EIGEN_FFTW_DEFAULT


System::System(Wavefunction wf, Potential pot) {
  matrixContraction = { Eigen::IndexPair<int>(1, 0) };
  rows = {0};
  columns = {1};
  addWavefunction(wf);
  addPotential(pot, 0, 0);
  timeStep = 0.0;
  threshold = 1e-20;
  nChannel = 1;
}

System::~System() {
//  std::cout << "System deleted" << std::endl;
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
  nChannel += 1;
}

void System::addZeroWavefunction(const double &ReducedMass, const double &Epsilon) {
  Wavefunction wf(wavefunctions[0].grid, ReducedMass, Epsilon);
  wf.initZero();
  wavefunctions.push_back(wf);
  times.push_back(doubleVec());
  energies.push_back(doubleVec());
  norms.push_back(doubleVec());
  averages.push_back(doubleVec());
  nChannel += 1;
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
//    cdMatrix potL(wavefunctions[0].grid.nPoint, 1);
//    if (potLeft[l] == potRight[l]) {
//      potL = (potentials[l].V + wavefunctions[potLeft[l]].epsilon).matrix();
//    }
//    else {
//      potL = potentials[l].V.matrix();
//    }
    cdMatrix potL = potentials[l].V.matrix();
    (potentialTensor.chip(potLeft[l],0)).chip(potRight[l],0) = Matrix_to_Tensor(potL, wavefunctions[0].grid.nPoint);
  }
}

void System::addGaussianPotential(const double &xCentre, const cd &height, const cd &sigma, const int &j, const int &k){
  Potential pot(wavefunctions[0].grid);
  pot.initZero();
  pot.addGaussian(xCentre, height, sigma);
  addPotential(pot, j, k);
}

void System::addConstantPotential(const cd &c, const double &xmin, const double &xmax, const int &j, const int &k){
  Potential pot(wavefunctions[0].grid);
  pot.initZero();
  pot.addConstant(c, xmin, xmax);
  addPotential(pot, j, k);
}

void System::addParabolicPotential(const double &xCentre, const cd &c, const int &j, const int &k){
  Potential pot(wavefunctions[0].grid);
  pot.initZero();
  pot.addParabolic(xCentre, c);
  addPotential(pot, j, k);
}

void System::evolveStep(int index, double timeStep, int maxOrder){
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

void System::evolveAllStep(double timeStep, int maxOrder){
  // Evolve all wavefunctions one step
  for (int j = 0; j < wavefunctions.size(); ++j) {
    evolveStep(j, timeStep, maxOrder);
  }
}

void System::evolveAll(int nSteps, double timeStep, int maxOrder) {
  // Evolve all wavefunctions by nSteps
  for (int j = 0; j < nSteps; ++j) {
    evolveAllStep(timeStep, maxOrder);
  }
}

void System::initCC(double tStep) {
  /// Setting up tensor to apply to wavefunctions
  // Assuming all based on the same grid, with same size
  // Initialise with the timestep that will be used for evolution
  timeStep = tStep;
  nChannel = unsigned(wavefunctions.size());
  psiTensor = cdVectorTensor(nChannel, wavefunctions[0].grid.nPoint);
  psiTensor.setZero();
  for (int j = 0; j < nChannel; ++j){
    cdMatrix psiMatrix = wavefunctions[j].psi.matrix();
    psiTensor.chip(j,0) = Matrix_to_Tensor(psiMatrix, wavefunctions[j].grid.nPoint);
  }
  U = cdMatrixTensor(nChannel, nChannel, wavefunctions[0].grid.nPoint);
  U.setZero();
  Udagger = cdMatrixTensor(nChannel, nChannel, wavefunctions[0].grid.nPoint);
  Udagger.setZero();
  cdVectorTensor D = cdVectorTensor(nChannel, wavefunctions[0].grid.nPoint);
  D.setZero();
  // Diagonalisation for finding U, Udagger and expD
  for (int j = 0; j < wavefunctions[0].grid.nPoint; ++j){
    cdVectorTensor potChip = potentialTensor.chip(j,2);
    cdMatrix potMat = Tensor_to_Matrix(potChip, nChannel, nChannel);
    ComplexEigenSolver<MatrixXcd> ces;
    // Threshold values
    for (int k = 0; k < nChannel; ++k) {
      for (int l = 0; l < nChannel; ++l) {
        if(potMat.real()(k,l) < threshold){
          potMat.real()(k,l) = 0.0;
        }
        if(potMat.imag()(k,l) < threshold){
          potMat.imag()(k,l) = 0.0;
        }
      }
    }
    ces.compute(potMat);
    U.chip(j,2) = Matrix_to_Tensor(ces.eigenvectors(), nChannel, nChannel) ;// U
    D.chip(j,1) = Vector_to_Tensor(ces.eigenvalues(), nChannel) ;// D
    cdMatrix Uinv = ces.eigenvectors().inverse();
    Udagger.chip(j,2) = Matrix_to_Tensor(Uinv, nChannel, nChannel) ;// Udagger
  }
  expD = ((-i*timeStep*0.5/HBARC)*D).exp();
  // Fourier Transform of Laplacian is momentum operator.
  for (int j = 0; j < nChannel; ++j){
    expP.emplace_back(exp((-i*timeStep*HBARC/(2.0*wavefunctions[j].reducedMass))*square(wavefunctions[j].grid.k)));
  }
  potentialOperator = cdMatrixTensor(nChannel, nChannel, wavefunctions[0].grid.nPoint);
  cdMatrixTensor UexpD = cdMatrixTensor(nChannel, nChannel, wavefunctions[0].grid.nPoint);
  for (int row = 0; row < nChannel; ++row) {
    UexpD.chip(row, 0) = (U.chip(row, 0))*expD;
  }
  for (int k = 0; k < nChannel; ++k) {
    for (int j = 0; j < nChannel; ++j) {
      potentialOperator.chip(k,0).chip(j,0) = (UexpD.chip(k,0)*Udagger.chip(j,1)).sum(rows);
    }
  }
  for (int k = 0; k < nChannel; ++k) {
    potentialOperator.chip(k,0) = (std::exp(-i*timeStep*(wavefunctions[k].epsilon)*0.25/HBARC))*potentialOperator.chip(k,0);
  }
  for (int k = 0; k < nChannel; ++k) {
    potentialOperator.chip(k,1) = (std::exp(-i*timeStep*(wavefunctions[k].epsilon)*0.25/HBARC))*potentialOperator.chip(k,1);
  }
  // Save the initial wavefunctions for Transmission computations
  for (auto wf: wavefunctions){
    wf.computePsiK();
    initialPsiKs.push_back(wf.psiK);
  }
}

void System::evolveCCStep(){
  // Apply half potential step
  cdVectorTensor step = cdVectorTensor(nChannel ,wavefunctions[0].grid.nPoint);
  for (int j = 0; j < nChannel; ++j) {
    step.chip(j,0) = (potentialOperator.chip(j,0)*psiTensor).sum(rows);
  }
  step.chip(0,1).setZero();
  step.chip(1,1).setZero();
  step.chip(wavefunctions[0].grid.nStep-1,1).setZero();
  step.chip(wavefunctions[0].grid.nStep,1).setZero();
// Apply kinetic step
  for (int k = 0; k < nChannel; ++k) {
    Eigen::Tensor<cd, 1> step1Chip = step.chip(k,0);
    cdVector psi = Tensor_to_Vector(step1Chip, wavefunctions[0].grid.nPoint);
    // FFT
    cdVector psi_input = (psi.array()*((-1.0*i*wavefunctions[0].grid.kMin*wavefunctions[0].grid.x).exp())).matrix();
    Eigen::FFT<double> fft;
    cdVector psi_output;
    psi_output.setZero(wavefunctions[0].grid.nPoint);
    fft.fwd(psi_output, psi_input);
    cdArray psiK = (psi_output.array())*(exp(-1.0*i*wavefunctions[0].grid.xMin*wavefunctions[0].grid.k))*wavefunctions[0].grid.xStep/(sqrt(2.0 * M_PI));
    // Apply momentum operator
    psiK *= expP[k];
    // IFFT
    psi_input = (psiK*((i*wavefunctions[0].grid.xMin*wavefunctions[0].grid.k).exp())*(sqrt(2.0* M_PI) / wavefunctions[0].grid.xStep)).matrix();
    psi_output.setZero(wavefunctions[0].grid.nPoint);
    fft.inv(psi_output, psi_input);
    //
    psi = ((psi_output.array())*(exp(i*wavefunctions[0].grid.kMin*wavefunctions[0].grid.x))).matrix();
    step.chip(k,0) = Vector_to_Tensor(psi, wavefunctions[0].grid.nPoint);
  }
  step.chip(0,1).setZero();
  step.chip(1,1).setZero();
  step.chip(wavefunctions[0].grid.nStep-1,1).setZero();
  step.chip(wavefunctions[0].grid.nStep,1).setZero();
  // Apply other half potential step
  for (int j = 0; j < nChannel; ++j) {
    psiTensor.chip(j,0) = (potentialOperator.chip(j,0)*step).sum(rows);
  }
  step.chip(0,1).setZero();
  step.chip(1,1).setZero();
  step.chip(wavefunctions[0].grid.nStep-1,1).setZero();
  step.chip(wavefunctions[0].grid.nStep,1).setZero();
}

void System::evolveCC(int nSteps) {
  cout << "Starting new system" << endl;
  boost::timer::progress_display show_progress(nSteps);
  for (int j = 0; j < nSteps; ++j) {
    evolveCCStep();
    ++show_progress;
  }
}

void System::evolveCC(int nSteps, std::vector<double> energies, std::vector<std::vector<double>> &data) {
  cout << "Starting new system" << endl;
  std::vector<double> timeVector;
  data.push_back(timeVector);
  for (int l = 0; l < energies.size(); ++l) {
    std::vector<double> totVec;
    data.push_back(totVec);
    for (int k = 0; k < nChannel; ++k) {
      std::vector<double> wfVector;
      data.push_back(wfVector);
    }
  }
  boost::timer::progress_display show_progress(nSteps);
  for (int j = 0; j < nSteps; ++j) {
    evolveCCStep();
    ++show_progress;
//    if ((j%int(0.5+5/timeStep)==0 and j*timeStep < 500) or (j%int(0.5+100/timeStep)==0))
    if (j%int(0.5+1/timeStep)==0){
      updateFromCC();
      for (int k = 0; k < nChannel; ++k) {
        times[k].push_back(j*timeStep);
        norms[k].push_back(wavefunctions[k].getNorm());
        averages[k].push_back(wavefunctions[k].getAvgX());
      }
      data[0].push_back(j*timeStep);
      // HELP
      for (int l = 0; l < energies.size(); ++l) {
        int dataIndex = (nChannel+1)*l+1;
        double T = 0.0;
        for (auto wf: wavefunctions) {
          double kPrime = std::sqrt(2 * wf.reducedMass * (energies[l] - wf.epsilon)) / HBARC;
          int before;
          int after;
          for (int k = int(1 + (wf.grid.nPoint) / 2); k < wf.grid.nPoint; ++k) {
            if (wf.grid.k(k) == kPrime) {
              T += norm(wf.psiK(k));
              break;
            } else if (wf.grid.k(k) > kPrime) {
              before = k - 1;
              after = k;
              T += norm(wf.psiK(before)) + (kPrime - wf.grid.k(before))
                  * ((norm(wf.psiK(after)) - norm(wf.psiK(before)))
                      / (wf.grid.k(after) - wf.grid.k(before))); // Add interpolated value
              break;
            }
          }
        }
        double kPrime = std::sqrt(2 * wavefunctions[0].reducedMass * (energies[l])) / HBARC;
        int before;
        int after;
        for (int k = int(1 + (wavefunctions[0].grid.nPoint) / 2); k < wavefunctions[0].grid.nPoint; ++k) {
          if (wavefunctions[0].grid.k(k) == kPrime) {
            T = T/norm(initialPsiKs[0](k));
            break;
          } else if (wavefunctions[0].grid.k(k) > kPrime) {
            before = k - 1;
            after = k;
            double orig = norm(initialPsiKs[0](before)) + (kPrime - wavefunctions[0].grid.k(before))
                * ((norm(initialPsiKs[0](after)) - norm(initialPsiKs[0](before)))
                    / (wavefunctions[0].grid.k(after) - wavefunctions[0].grid.k(before))); // Add interpolated value
            T = T/orig;
            break;
          }
        }
        // Find reflections
        std::vector<double> Rs(nChannel);
        for (int n = 0; n < nChannel; ++n) {
          double R = 0.0;
          kPrime = std::sqrt(2 * wavefunctions[n].reducedMass * (energies[l] - wavefunctions[n].epsilon)) / HBARC;
          for (int k = int(1 + (wavefunctions[n].grid.nPoint) / 2); k < wavefunctions[n].grid.nPoint; ++k) {
            if (wavefunctions[n].grid.k(k) == kPrime) {
              int minus = wavefunctions[n].grid.nPoint - k;
              R += norm(wavefunctions[n].psiK(minus));
              break;
            } else if (wavefunctions[n].grid.k(k) > kPrime) {
              before = wavefunctions[n].grid.nPoint - (k-1);
              after = wavefunctions[n].grid.nPoint  - k;
              R += norm(wavefunctions[n].psiK(before)) + (-kPrime - wavefunctions[n].grid.k(before))
                  * ((norm(wavefunctions[n].psiK(after)) - norm(wavefunctions[n].psiK(before)))
                      / (wavefunctions[n].grid.k(after) - wavefunctions[n].grid.k(before))); // Add interpolated value
              break;
            }
          }
          double kPrime = std::sqrt(2 * wavefunctions[0].reducedMass * (energies[l])) / HBARC;
          int before;
          int after;
          for (int k = int(1 + (wavefunctions[0].grid.nPoint) / 2); k < wavefunctions[0].grid.nPoint; ++k) {
            if (wavefunctions[0].grid.k(k) == kPrime) {
              R = R/norm(initialPsiKs[0](k));
              Rs[n] = R;
              break;
            } else if (wavefunctions[0].grid.k(k) > kPrime) {
              before = k - 1;
              after = k;
              double orig = norm(initialPsiKs[0](before)) + (kPrime - wavefunctions[0].grid.k(before))
                  * ((norm(initialPsiKs[0](after)) - norm(initialPsiKs[0](before)))
                      / (wavefunctions[0].grid.k(after) - wavefunctions[0].grid.k(before))); // Add interpolated value
              R = R/orig;
              Rs[n] = R;
              break;
            }
          }
        }
        double Rtot = std::accumulate(Rs.begin(), Rs.end(), double(0.0));
        data[dataIndex].push_back(T+Rtot);
        for (int m = 0; m < Rs.size(); ++m){
          data[dataIndex+m+1].push_back(Rs[m]);
        }
      }
    }
  }
  updateFromCC();
}

void System::updateFromCC(){
  for (int k = 0; k < nChannel; ++k) {
    Eigen::Tensor<cd, 1> psiTensorChip = psiTensor.chip(k,0);
    cdVector psi = Tensor_to_Vector(psiTensorChip, wavefunctions[0].grid.nPoint);
    wavefunctions[k].psi = psi.array();
    wavefunctions[k].computePsiK();
  }
}

void System::updateK(){
  for (int j = 0; j < nChannel; ++j) {
    wavefunctions[j].computePsiK();
  }
}

void System::log(double time){
  if(timeStep != 0.0){
    updateFromCC();
  }
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

dArray System::getTransmission(){
  /// TODO: Consider using Splinter for interpolation (easy splines to sample, but requires compilation)
  dArray T;
  T.setZero(wavefunctions[0].grid.nPoint);
  for (auto wf: wavefunctions){
    if (wf.epsilon==0.0){
      for (int m = int(1+(wf.grid.nPoint)/2); m < wf.grid.nPoint; ++m) {
        T(m) += wf.psiK.abs2()(m);
      }
    }
    else {
      for (int m = int(1+(wf.grid.nPoint)/2); m < wf.grid.nPoint; ++m) {
        double E = wavefunctions[0].E(m);
        if (E - wf.epsilon >= 0.0){
          double kPrime = std::sqrt(2*wf.reducedMass*(E-wf.epsilon))/HBARC;
          int before;
          int after;
          for (int k = int(1+(wf.grid.nPoint)/2); k < wf.grid.nPoint; ++k) {
            if (wf.grid.k(k) == kPrime){
              T(m) += wf.psiK.abs2()(k);
              break;
            }
            else if (wf.grid.k(k) > kPrime){
              before = k-1;
              after = k;
              T(m) += wf.psiK.abs2()(before) + (kPrime - wf.grid.k(before))*((wf.psiK.abs2()(after)-wf.psiK.abs2()(before))/(wf.grid.k(after)-wf.grid.k(before))); // Add interpolated value
              break;
            }
          }
        }
      }
    }
  }
  T = T/initialPsiKs[0].abs2();
  return T;
}

dArray System::getReflection(int index){
  /// TODO: Consider using Splinter for interpolation (easy splines to sample, but requires compilation)
  dArray R;
  R.setZero(wavefunctions[0].grid.nPoint);
  if (wavefunctions[index].epsilon==0.0){
    for (int m = int(1+(wavefunctions[index].grid.nPoint)/2); m < wavefunctions[index].grid.nPoint; ++m) {
      int minus = wavefunctions[index].grid.nPoint - m;
      R(m) += wavefunctions[index].psiK.abs2()(minus);
    }
  }
  else {
    for (int m = int(1+(wavefunctions[index].grid.nPoint)/2); m < wavefunctions[index].grid.nPoint; ++m) {
      double E = wavefunctions[0].E(m);
      if (E - wavefunctions[index].epsilon >= 0.0){
        double kPrime = std::sqrt(2*wavefunctions[index].reducedMass*(E-wavefunctions[index].epsilon))/HBARC;
        int before;
        int after;
        for (int k = int(1+(wavefunctions[index].grid.nPoint)/2); k < wavefunctions[index].grid.nPoint; ++k) {
          if (wavefunctions[index].grid.k(k) == kPrime){
            int minus = wavefunctions[index].grid.nPoint - k;
            R(m) += wavefunctions[index].psiK.abs2()(minus);
            break;
          }
          else if (wavefunctions[index].grid.k(k) > kPrime){
            before = wavefunctions[index].grid.nPoint - (k-1);
            after = wavefunctions[index].grid.nPoint  - k;
            R(m) += wavefunctions[index].psiK.abs2()(before) + (-kPrime - wavefunctions[index].grid.k(before))*((wavefunctions[index].psiK.abs2()(after)-wavefunctions[index].psiK.abs2()(before))/(wavefunctions[index].grid.k(after)-wavefunctions[index].grid.k(before))); // Add interpolated value
            break;
          }
        }
      }
    }
  }
  R = R/initialPsiKs[0].abs2();
  return R;
}
