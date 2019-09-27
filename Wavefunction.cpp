#include <iostream>
#include <vector>
#include <complex>
#include "eigen/Eigen/Core"
#include "eigen/unsupported/Eigen/FFT"
#define EIGEN_FFTW_DEFAULT
#include "Wavefunction.h"

Wavefunction::Wavefunction(const Grid &object, const double &ReducedMass, const double &Epsilon) : grid(1, 0.0, 1.0, 1.0) {
  grid = object;
  reducedMass = ReducedMass * AMU;
  psiK.resize(grid.nPoint, Eigen::NoChange);
  psi.resize(grid.nPoint, Eigen::NoChange);
  epsilon = Epsilon;
}

Wavefunction::Wavefunction(const Grid &object, const double &ReducedMass) : grid(1, 0.0, 1.0, 1.0) {
  grid = object;
  reducedMass = ReducedMass * AMU;
  psiK.resize(grid.nPoint, Eigen::NoChange);
  psi.resize(grid.nPoint, Eigen::NoChange);
  epsilon = 0.0;
  E = ((HBARC*grid.k).abs2())/(2.0*reducedMass);
}

Wavefunction::~Wavefunction() {
//  std::cout << "Wavefunction deleted" << std::endl;
}

void Wavefunction::test() {
  std::cout << "Test Wavefunction" << std::endl;
  std::cout << "Psi is: " << psi << std::endl;
}


void Wavefunction::normalise() {
  psi = (1.0/sqrt(getNorm()))*psi;
}

void Wavefunction::initZero() {
  psi.setZero(grid.nPoint);
}

void Wavefunction::initGaussian(const double &mean, const double &sigma) {
  psi = exp(-1.0*(square(grid.x - mean))/(2.0*sigma*sigma));
  zeroEdges();
  normalise();
}

void Wavefunction::initAsymmGaussian(const double &mean, const double &sigma) {
  psi = exp(-1.0*(square(grid.x - mean - sigma))/(2.0*sigma*sigma)) - exp(-1.0*(square(grid.x - mean + sigma))/(2.0*sigma*sigma));
  zeroEdges();
  normalise();
}

void Wavefunction::zeroEdges() {
  psi(0) = cd(0.0, 0.0);
  psi(grid.nStep) = cd(0.0, 0.0);
}

void Wavefunction::initSine(const double &N) {
  psi = (N*M_PI*(grid.x-grid.xMin)/(grid.xMax-grid.xMin)).sin();
  normalise();
}

void Wavefunction::initConstant() {
  psi.setConstant(1.0);
  zeroEdges();
  normalise();
}

void Wavefunction::initConstantInRegion(const double &xmin, const double &xmax) {
  psi.setZero(grid.nPoint);
  for (int j = 0; j < grid.nPoint; ++j) {
    if (grid.x(j) > xmin and grid.x(j) < xmax) {
      psi(j) = 1.0;
    }
  }
  zeroEdges();
  normalise();
}

void Wavefunction::boostWaveNumber(const double &WN) {
  psi *= (exp(i*WN*grid.x));
  normalise();
}

void Wavefunction::boostEnergy(const double &energy) {
  double WN = sgn(energy) * sqrt(2 * reducedMass * std::abs(energy)) / HBARC;
  boostWaveNumber(WN);
}

/// Getters

double Wavefunction::getNorm() {
  return overlap((*this));
}

double Wavefunction::getNormInRegion(const double &xmin, const double &xmax) {
  doubleVec integrand;
  integrand.reserve(grid.nPoint);
  for (int j = 0; j < grid.nPoint; ++j) {
    if (grid.x(j) > xmin and grid.x(j) < xmax) {
      integrand.push_back(std::norm(psi(j)));
    }
  }
  assert(("No points in this range", !integrand.empty()));
  return vectorTrapezoidIntegrate(integrand, grid.xStep, int(integrand.size())-1);
}

dArray Wavefunction::getReal() {
  return psi.real();
}

dArray Wavefunction::getImag() {
  return psi.imag();
}

dArray Wavefunction::getAbs() {
  return psi.abs();
}

dArray Wavefunction::getAbsSq() {
  return psi.abs2();
}

dArray Wavefunction::getKReal() {
  return psiK.real();
}

dArray Wavefunction::getKImag() {
  return psiK.imag();
}

dArray Wavefunction::getKAbs() {
  return psiK.abs();
}

dArray Wavefunction::getKAbsSq() {
  return psiK.abs2();
}

double Wavefunction::getAvgX() {
  dArray integrand = (psi.abs2())*grid.x;
  return vectorTrapezoidIntegrate(integrand, grid.xStep, grid.nStep);
}

void Wavefunction::copy(const Wavefunction &wf) {
  grid = wf.grid;
  reducedMass = wf.reducedMass;
  psi = wf.psi;
  psiK = wf.psiK;
}

double Wavefunction::overlap(const Wavefunction &object) {
  dArray integrand = psi.abs2();
  return vectorTrapezoidIntegrate(integrand, grid.xStep, grid.nStep);
}

/// TODO: Fix introduced error of ~1e-16 in computations

void Wavefunction::computePsi() {
  cdVector psi_input = (psiK*((i*grid.xMin*grid.k).exp())*(sqrt(2.0* M_PI) / grid.xStep)).matrix();
  // FFT
  Eigen::FFT<double> fft;
  cdVector psi_output;
  psi_output.setZero(grid.nPoint);
  fft.inv(psi_output, psi_input);
  //
  psi = (psi_output.array())*(exp(i*grid.kMin*grid.x));
}

void Wavefunction::computePsiK(){
  cdVector psi_input = (psi*((-1.0*i*grid.kMin*grid.x).exp())).matrix();
  // FFT
  Eigen::FFT<double> fft;
  cdVector psi_output;
  psi_output.setZero(grid.nPoint);
  fft.fwd(psi_output, psi_input);
  //
  psiK = (psi_output.array())*(exp(-1.0*i*grid.xMin*grid.k))*grid.xStep/(sqrt(2.0 * M_PI));
}


//void Wavefunction::computePsiK() {
//  // Compute input for FFT
//  cdVector step1 = (psi*(exp(-1.0*i*grid.kMin*grid.x))).matrix();
//  complexVec psi_input;
//  psi_input.resize(step1.size());
//  cdVector::Map(&psi_input[0], step1.size()) = step1;
//  doubleVec dvec = fourierComplexToDouble(psi_input);  // Need input as a double array to Fourier Transform
//
//  // FFT
//  gsl_fft_complex_wavetable *wavetable;
//  gsl_fft_complex_workspace *workspace;
//  wavetable = gsl_fft_complex_wavetable_alloc(grid.nPoint);
//  workspace = gsl_fft_complex_workspace_alloc(grid.nPoint);
//  int res = gsl_fft_complex_forward(dvec.data(), 1, grid.nPoint, wavetable, workspace);
//  if (res != 0) {
//    std::cout << "FFT Failed" << std::endl;
//  }
//  gsl_fft_complex_wavetable_free(wavetable);
//  gsl_fft_complex_workspace_free(workspace);
//
//  std::rotate(dvec.begin(), dvec.begin() + int(dvec.size()/2), dvec.end());
//  psi_input = fourierDoubleToComplex(dvec);
//  cdArray psi_output = Map<cdArray>(psi_input.data(), psi_input.size());
//  // Convert back
//  psiK = (psi_output)*((-1.0*i*grid.xMin*grid.k).exp())*grid.xStep/(sqrt(2.0 * M_PI));
//}

//void Wavefunction::computePsi() {
//  // Compute input for FFT
//  cdVector step1 = (psiK*((i*grid.xMin*grid.k).exp())*(sqrt(2.0* M_PI) / grid.xStep)).matrix();
//  complexVec psi_input;
//  psi_input.resize(step1.size());
//  cdVector::Map(&psi_input[0], step1.size()) = step1;
//  std::rotate(psi_input.begin(), psi_input.begin()+int((psi_input.size()+1)/2), psi_input.end());  // Rotate to match ordering of FFT
//  doubleVec dvec = fourierComplexToDouble(psi_input);  // Need input as a double array to Fourier Transform
//
//  // IFFT
//  gsl_fft_complex_wavetable *wavetable;
//  gsl_fft_complex_workspace *workspace;
//  wavetable = gsl_fft_complex_wavetable_alloc(grid.nPoint);
//  workspace = gsl_fft_complex_workspace_alloc(grid.nPoint);
//  int res = gsl_fft_complex_inverse(dvec.data(), 1, grid.nPoint, wavetable, workspace);
//  if (res != 0) {// Check IFFT worked
//    std::cout << "FFT Failed" << std::endl;
//  }
//  gsl_fft_complex_wavetable_free(wavetable);
//  gsl_fft_complex_workspace_free(workspace);
//
//  psi_input = fourierDoubleToComplex(dvec);
//  cdArray psi_output = Map<cdArray>(psi_input.data(), psi_input.size());
//  // Convert back
//  psi = (psi_output)*((i*grid.kMin*grid.x).exp());
//}