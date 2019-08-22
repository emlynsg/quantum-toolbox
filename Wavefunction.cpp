#include <iostream>
#include <vector>
#include <complex>
#include "eigen/Eigen/Core"
#include "eigen/unsupported/Eigen/FFT"

#include "Wavefunction.h"

Wavefunction::Wavefunction(const Grid &object, const double &ReducedMass) : grid(1, 0.0, 1.0, 1.0) {
  grid = object;
  reducedMass = ReducedMass * AMU;
  psiK.resize(grid.nPoint, Eigen::NoChange);
  psi.resize(grid.nPoint, Eigen::NoChange);
}

Wavefunction::~Wavefunction() {
  std::cout << "Wavefunction deleted" << std::endl;
}

void Wavefunction::test() {
  std::cout << "Test Wavefunction" << std::endl;
  std::cout << "Average value: " << getAvgX() << std::endl;
}


void Wavefunction::normalise() {
  psi = (1.0/sqrt(getNorm()))*psi;
}


/// Check ordering on PsiK and Psi outputs
/// https://www.gnu.org/software/gsl/doc/html/fft.html#overview-of-complex-data-ffts
/// Note: Definitely wrong as ifft then fft doesn't give the original values


void Wavefunction::computePsiK() {
  // Compute input for FFT
  complexVec psi_input = vectorMultiply(psi, vectorExp(vectorScale(grid.x, -1. * grid.kMin * i)));
  doubleVec dvec = fourierComplexToDouble(psi_input);  // Need input as a double array to Fourier Transform


  // FFT
  gsl_fft_complex_wavetable *wavetable;
  gsl_fft_complex_workspace *workspace;
  wavetable = gsl_fft_complex_wavetable_alloc(grid.nPoint);
  workspace = gsl_fft_complex_workspace_alloc(grid.nPoint);
  int res = gsl_fft_complex_forward(dvec.data(), 1, grid.nPoint, wavetable, workspace);
  if (res != 0) {
    std::cout << "FFT Failed" << std::endl;
  }
  gsl_fft_complex_wavetable_free(wavetable);
  gsl_fft_complex_workspace_free(workspace);

  std::rotate(dvec.begin(), dvec.begin() + int(dvec.size()/2), dvec.end());
  // Convert back
  psiK = vectorScale(vectorMultiply(fourierDoubleToComplex(dvec), vectorExp(vectorScale(grid.k, -1. * grid.xMin * i))),
                     grid.xStep / (sqrt(2. * M_PI)));
}

void Wavefunction::computePsi() {
  // Compute input for FFT
  complexVec psi_input =
      vectorScale(vectorMultiply(psiK, vectorExp(vectorScale(grid.k, grid.xMin * i))), (sqrt(2. * M_PI) / grid.xStep));
  std::rotate(psi_input.begin(), psi_input.begin()+int((psi_input.size()+1)/2), psi_input.end());  // Rotate to match ordering of FFT
  doubleVec dvec = fourierComplexToDouble(psi_input);  // Need input as a double array to Fourier Transform

  // IFFT
  gsl_fft_complex_wavetable *wavetable;
  gsl_fft_complex_workspace *workspace;
  wavetable = gsl_fft_complex_wavetable_alloc(grid.nPoint);
  workspace = gsl_fft_complex_workspace_alloc(grid.nPoint);
  int res = gsl_fft_complex_inverse(dvec.data(), 1, grid.nPoint, wavetable, workspace);
  if (res != 0) {// Check IFFT worked
    std::cout << "FFT Failed" << std::endl;
  }
  gsl_fft_complex_wavetable_free(wavetable);
  gsl_fft_complex_workspace_free(workspace);

  // Convert back
  psi = vectorMultiply(fourierDoubleToComplex(dvec), vectorExp(vectorScale(grid.x, 1. * grid.kMin * i)));
}

void Wavefunction::initZero() {
  psi.setZero(grid.nPoint);
}

void Wavefunction::initGaussian(const double &mean, const double &sigma) {
  psi = (((grid.x - mean).square()).exp())/(2*sigma*sigma);
  zeroEdges();
  normalise();
}

void Wavefunction::initAsymmGaussian(const double &mean, const double &sigma) {
  psi = ((((grid.x-mean-sigma).square()).exp()) - (((grid.x-mean+sigma).square()).exp()))/(2*sigma*sigma);
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

void Wavefunction::boostWaveNumber(const double &WN) {
  psi *= (i*WN*grid.x).exp();
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
    if (grid.x[j] >= xmin and grid.x[j] <= xmax) {
      integrand.push_back(std::norm(psi[j]));
    }
  }
  assert(("No points in this range", !integrand.empty()));
  return vectorTrapezoidIntegrate(integrand, grid.xStep, int(integrand.size()));
}

dVec Wavefunction::getReal() {
  return psi.real();
}

dVec Wavefunction::getImag() {
  return psi.imag();
}

dVec Wavefunction::getAbs() {
  return psi.abs();
}

dVec Wavefunction::getAbsSq() {
  return psi.abs2();
}

dVec Wavefunction::getKReal() {
  return psiK.real();
}

dVec Wavefunction::getKImag() {
  return psiK.imag();
}

dVec Wavefunction::getKAbs() {
  return psiK.abs();
}

dVec Wavefunction::getKAbsSq() {
  return psiK.abs2();
}

double Wavefunction::getAvgX() {
  dVec integrand = (psi.abs2())*grid.x;
  return vectorTrapezoidIntegrate(integrand, grid.xStep, grid.nPoint);
}

void Wavefunction::copy(const Wavefunction &wf) {
  grid = wf.grid;
  reducedMass = wf.reducedMass;
  psi = wf.psi;
  psiK = wf.psiK;
}

double Wavefunction::overlap(const Wavefunction &object) {
  dVec integrand = psi.abs2();
  return vectorTrapezoidIntegrate(integrand, grid.xStep, grid.nPoint);
}