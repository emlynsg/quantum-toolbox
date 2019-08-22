#include "Extras.h"

/// Variables

cd i = cd(0.0, 1.0);
double HBARC = 197.3; // MeV fm
double AMU = 931.5;   // MeV/c^2
double ESQ = 1.44;    // Electron charge in MeV fm

/// Distribution functions

double sine(const double &x, const double &N, const double &xmin, const double &xmax) {
  return sin(N * M_PI * (x - xmin) / (xmax - xmin));
}

double gaussian(const double &x, const double &X0, const double &Sigma) {
  return exp(-pow(x - X0, 2.0) / (2 * Sigma * Sigma));
}

double asymmGaussian(const double &x, const double &X0, const double &Sigma) {
  return 1.0
      * (exp(-pow(x - X0 - Sigma, 2.0) / (2 * Sigma * Sigma)) - exp(-pow(x - X0 + Sigma, 2.0) / (2 * Sigma * Sigma)));
}


/// Vector functions

doubleVec vectorComplexToDouble(const complexVec &a) {
  doubleVec b(a.size());
  std::transform(a.begin(), a.end(), b.begin(), [](complex elt) { return elt.real(); });
  return b;
}

complexVec vectorDoubleToComplex(const doubleVec &a) {
  complexVec b(a.size());
  std::transform(a.begin(), a.end(), std::back_inserter(b),
                 [](double r) { return std::complex<double>(r, 0.0); });
  return b;
}

complexVec vectorMultiply(const complexVec &a, const complexVec &b) {
  assert(("Vector lengths don't match", a.size() == b.size()));
  complexVec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::multiplies<>());
  return c;
}

doubleVec vectorMultiply(const doubleVec &a, const doubleVec &b) {
  assert(("Vector lengths don't match", a.size() == b.size()));
  doubleVec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::multiplies<>());
  return c;
}

complexVec vectorExp(const complexVec &a) {
  complexVec b(a.size());
  std::transform(a.begin(), a.end(), b.begin(), expVec());
  return b;
}

complexVec vectorScale(const complexVec &a, const complex &b) {
  complexVec c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto &elt) { return elt * b; });
  return c;
}

complexVec vectorScale(const complexVec &a, const double &b) {
  complexVec c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto &elt) { return elt * b; });
  return c;
}

doubleVec vectorScale(const doubleVec &a, const double &b) {
  doubleVec c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto &elt) { return elt * b; });
  return c;
}

complexVec vectorScale(const doubleVec &a, const complex &b) {
  complexVec c(a.size());
  std::transform(a.begin(), a.end(), c.begin(), [b](auto &elt) { return elt * b; });
  return c;
}

complexVec vectorAdd(const complexVec &a, const complexVec &b) {
  complexVec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<>());
  return c;
}

complexVec vectorAdd(const complexVec &a, const doubleVec &b) {
  complexVec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<>());
  return c;
}

complexVec vectorAdd(const doubleVec &a, const complexVec &b) {
  complexVec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<>());
  return c;
}

complexVec vectorSubtract(const complexVec &a, const complexVec &b) {
  complexVec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::minus<>());
  return c;
}

complexVec vectorSubtract(const complexVec &a, const doubleVec &b) {
  complexVec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::minus<>());
  return c;
}

complexVec vectorSubtract(const doubleVec &a, const complexVec &b) {
  complexVec c(a.size());
  std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::minus<>());
  return c;
}

doubleVec fourierComplexToDouble(const complexVec &cvector) {
  doubleVec dvector;
  dvector.reserve((2 * cvector.size()));
  for (auto j : cvector) {
    dvector.push_back(j.real());
    dvector.push_back(j.imag());
  }
  return dvector;
}

dVec fourierComplexToDouble(cVec &cvector) {
  dVec dvector;
  dvector.resize(2*cvector.size());
  for (int k = 0; k < cvector.size(); ++k) {
    dvector[2*k] = (cvector[k].real());
    dvector[2*k+1] = (cvector[k].imag());
  }
  return dvector;
}

complexVec fourierDoubleToComplex(const doubleVec &dvector) {
  complexVec cvector;
  cvector.reserve((dvector.size()) / 2);
  int range = 0;
  while (range < dvector.size()) {
    complex c = complex(dvector[range], dvector[range + 1]);
    cvector.push_back(c);
    range = range + 2;
  }
  return cvector;
}

cVec fourierDoubleToComplex(dVec &dvector) {
  cVec cvector;
  cvector.resize(dvector.size()/2);
  for (int k = 0; k < 2*cvector.size(); k=k+2) {
    cvector[k] = cd(dvector[k], dvector[k+1]);
  }
  return cvector;
}

double vectorTrapezoidIntegrate(const doubleVec &vect, const double &h, const int &n) {
  assert(("Integration requires a minimum of 9 points", n > 9));
  return (h / 2.0) * (vect[0] + 2.0 * std::accumulate(vect.begin() + 1, vect.begin() + (vect.size() - 1), 0.0)
      + vect[n]);
}

double vectorTrapezoidIntegrate(dVec &vect, const double &h, const int &n) {
  assert(("Integration requires a minimum of 9 points", n > 9));
  return (h / 2.0) * (vect[0] + 2.0 * vect.segment(1,n).sum() + vect[n]);
}

/// Simpson Rule (from Wikipedia, not sure of reference)
/// TODO: Fix integration approach

double vectorSimpsonIntegrate(const doubleVec &vect, const double &h, const int &n) {
  assert(("Integration requires a minimum of 9 points", n > 9));
  return (h / 48.0) * (17.0 * vect[0] + 59.0 * vect[1] + 43.0 * vect[2] + 49.0 * vect[3]
      + 48.0 * std::accumulate(vect.begin() + 4, vect.begin() + (vect.size() - 4), 0.0)
      + 49.0 * vect[n - 3] + 43.0 * vect[n - 2] + 59.0 * vect[n - 1] + 17.0 * vect[n]);
}

double vectorSimpsonIntegrate(dVec &vect, const double &h, const int &n) {
  assert(("Integration requires a minimum of 9 points", n > 9));
  return (h / 48.0) * (17.0 * vect[0] + 59.0 * vect[1] + 43.0 * vect[2] + 49.0 * vect[3]
      + 48.0 * vect.segment(4,n-3).sum()
      + 49.0 * vect[n - 3] + 43.0 * vect[n - 2] + 59.0 * vect[n - 1] + 17.0 * vect[n]);
}

