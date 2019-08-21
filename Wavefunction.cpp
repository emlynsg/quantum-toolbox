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
  //exp(-pow(x - X0, 2.0) / (2 * Sigma * Sigma));
  psi = (((grid.x - mean).square()).exp())/(2*sigma*sigma);
  zeroEdges();
  normalise();
}

void Wavefunction::initAsymmGaussian(const double &mean, const double &sigma) {
  1.0
      * (exp(-pow(x - X0 - Sigma, 2.0) / (2 * Sigma * Sigma)) - exp(-pow(x - X0 + Sigma, 2.0) / (2 * Sigma * Sigma)));
  psi = ((((grid.x-mean-sigma).square()).exp()) - (((grid.x-mean+sigma).square()).exp()))/(2*sigma*sigma);
  zeroEdges();
  normalise();
}

void Wavefunction::zeroEdges() {
  psi(0) = cd(0.0, 0.0);
  psi(grid.nStep) = cd(0.0, 0.0);
}

void Wavefunction::initSine(const double &N) {
  sin(N * M_PI * (x - xmin) / (xmax - xmin));
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
  for (int j = 0; j < grid.nPoint; ++j) {
    if (grid.x[j] >= xmin and grid.x[j] <= xmax) {
      integrand.push_back(std::norm(psi[j]));
    }
  }
  assert(("No points in this range", !integrand.empty()));
  return vectorTrapezoidIntegrate(integrand, grid.xStep, int(integrand.size()));
}

dVec Wavefunction::getReal() {
  doubleVec returnValue(grid.nPoint);
  std::transform(psi.begin(), psi.end(), returnValue.begin(), [](auto &elt) { return elt.real(); });
  return returnValue;
}

dVec Wavefunction::getImag() {
  doubleVec returnValue(grid.nPoint);
  std::transform(psi.begin(), psi.end(), returnValue.begin(), [](auto &elt) { return elt.imag(); });
  return returnValue;
}

dVec Wavefunction::getAbs() {
  doubleVec returnValue(grid.nPoint);
  std::transform(psi.begin(), psi.end(), returnValue.begin(), [](auto &elt) { return std::abs(elt); });
  return returnValue;
}

dVec Wavefunction::getAbsSq() {
  doubleVec returnValue(grid.nPoint);
  std::transform(psi.begin(), psi.end(), returnValue.begin(), [](auto &elt) { return pow(std::abs(elt), 2.0); });
  return returnValue;
}

dVec Wavefunction::getKReal() {
  doubleVec returnValue(grid.nPoint);
  std::transform(psiK.begin(), psiK.end(), returnValue.begin(), [](auto &elt) { return elt.real(); });
  return returnValue;
}

dVec Wavefunction::getKImag() {
  doubleVec returnValue(grid.nPoint);
  std::transform(psiK.begin(), psiK.end(), returnValue.begin(), [](auto &elt) { return elt.imag(); });
  return returnValue;
}

dVec Wavefunction::getKAbs() {
  doubleVec returnValue(grid.nPoint);
  std::transform(psiK.begin(), psiK.end(), returnValue.begin(), [](auto &elt) { return std::abs(elt); });
  return returnValue;
}

dVec Wavefunction::getKAbsSq() {
  doubleVec returnValue(grid.nPoint);
  std::transform(psiK.begin(), psiK.end(), returnValue.begin(), [](auto &elt) { return pow(std::abs(elt), 2.0); });
  return returnValue;
}

double Wavefunction::getAvgX() {
  doubleVec integrand(grid.nPoint);
  for (int j = 0; j < grid.nPoint; ++j) {
    integrand[j] = (std::abs(psi[j] * std::conj(psi[j])));
  }
  double returnValue = vectorTrapezoidIntegrate(vectorMultiply(integrand, grid.x), grid.xStep, grid.nPoint);
  return returnValue;
}

void Wavefunction::copy(const Wavefunction &wf) {
  grid = wf.grid;
  reducedMass = wf.reducedMass;
  psi.reserve(wf.grid.nPoint);
  psiK.reserve(wf.grid.nPoint);
  psi = wf.psi;
  psiK = wf.psiK;
}


/// TODO: Fix overlap by completing integration approach

double Wavefunction::overlap(const Wavefunction &object) {
  doubleVec integrand(grid.nPoint);
  for (int j = 0; j < grid.nPoint; ++j) {
    integrand[j] = (std::abs(psi[j] * std::conj(object.psi[j])));
  }
  double returnValue = vectorTrapezoidIntegrate(integrand, grid.xStep, grid.nPoint);
  return returnValue;
}



/// All GSL integration seems to need a function as input
/// TODO: Interpolate points, write this as a function, and then integrate

//double Wavefunction::overlap(const Wavefunction& object) {
//
//  struct spline_parameters{gsl_spline a; gsl_interp_accel b;};
//
//  double grid_array[grid.nPoint];
//  double integrand[grid.nPoint];
//  for(int j=0; j<grid.nPoint; ++j) {
//    grid_array[j] = grid.x[j];
//    integrand[j] = std::abs(psi[j] * std::conj(object.psi[j]));
//  }
//  gsl_interp_accel *acc = gsl_interp_accel_alloc();
//  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, grid.nPoint);
//  gsl_spline_init(spline, grid_array, integrand, grid.nPoint);
//
//  double integrand_value = [](double xi, void * parameters) {
//    struct spline_parameters *params = (struct spline_parameters *)parameters;
//    gsl_spline *splin = &(params->a);
//    gsl_interp_accel *ac = &(params->b);
//    double integrand_val = gsl_spline_eval(splin, xi, ac);
//    return integrand_val;
//  };
//
//  gsl_integration_workspace * w
//      = gsl_integration_workspace_alloc (1000);
//
//  double result, error;
//  gsl_function F;
//  struct spline_parameters params;
//  params.a = *spline;
//  params.b = *acc;
//  F.function = &integrand_value;
//  F.params = &params;
//
//  gsl_integration_qag(&F, grid.xMin, grid.xMax, 0, 1e-5
//                       , grid.nPoint, 6, w, &result, &error);
//
//  std::cout << "Result and error: " << result << ", " << error << std::endl;
//
//  gsl_spline_free (spline);
//  gsl_interp_accel_free (acc);
//
//  return result;
//}