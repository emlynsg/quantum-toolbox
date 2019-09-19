#include "Extras.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

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
  std::transform(a.begin(), a.end(), b.begin(), [](cd elt) { return elt.real(); });
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

complexVec vectorScale(const complexVec &a, const cd &b) {
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

complexVec vectorScale(const doubleVec &a, const cd &b) {
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

doubleVec fourierComplexToDouble(const complexVec &cdVector) {
  doubleVec dvector;
  dvector.reserve((2 * cdVector.size()));
  for (auto j : cdVector) {
    dvector.push_back(j.real());
    dvector.push_back(j.imag());
  }
  return dvector;
}

dArray fourierComplexToDouble(cdArray &cdvector) {
  dArray dvector;
  dvector.resize(2*cdvector.size());
  for (int k = 0; k < cdvector.size(); ++k) {
    dvector[2*k] = (cdvector[k].real());
    dvector[2*k+1] = (cdvector[k].imag());
  }
  return dvector;
}

complexVec fourierDoubleToComplex(const doubleVec &dvector) {
  complexVec cdvector;
  cdvector.reserve((dvector.size()) / 2);
  int range = 0;
  while (range < dvector.size()) {
    cd c = cd(dvector[range], dvector[range + 1]);
    cdvector.push_back(c);
    range = range + 2;
  }
  return cdvector;
}

cdArray fourierDoubleToComplex(dArray &dvector) {
  cdArray cdVector;
  cdVector.resize(dvector.size()/2);
  for (int k = 0; k < 2*cdVector.size(); k=k+2) {
    cdVector[k] = cd(dvector[k], dvector[k+1]);
  }
  return cdVector;
}

double vectorTrapezoidIntegrate(doubleVec &vect, const double &h, const int &n) {
  assert(("Integration requires a minimum of 9 points", n > 9));
  return (h / 2.0) * (vect[0] + 2.0 * std::accumulate(vect.begin() + 1, vect.begin() + (vect.size() - 1), 0.0)
      + vect[n]);
}

double vectorTrapezoidIntegrate(dArray &vect, const double &h, const int &n) {
  assert(("Integration requires a minimum of 9 points", n > 9));
  return (h / 2.0) * (vect[0] + 2.0 * vect.segment(1,n-1).sum() + vect[n]);
}

/// Simpson Rule (from Wikipedia, not sure of reference)
/// TODO: Fix integration approach

//double vectorSimpsonIntegrate(doubleVec &vect, double &h, int &n) {
//  assert(("Integration requires a minimum of 9 points", n > 9));
//  return (h / 48.0) * (17.0 * vect[0] + 59.0 * vect[1] + 43.0 * vect[2] + 49.0 * vect[3]
//      + 48.0 * std::accumulate(vect.begin() + 4, vect.begin() + (vect.size() - 4), 0.0)
//      + 49.0 * vect[n - 3] + 43.0 * vect[n - 2] + 59.0 * vect[n - 1] + 17.0 * vect[n]);
//}
//
//double vectorSimpsonIntegrate(dArray &vect, double &h, int &n) {
//  assert(("Integration requires a minimum of 9 points", n > 9));
//  return (h / 48.0) * (17.0 * vect[0] + 59.0 * vect[1] + 43.0 * vect[2] + 49.0 * vect[3]
//      + 48.0 * vect.segment(4,n-3).sum()
//      + 49.0 * vect[n - 3] + 43.0 * vect[n - 2] + 59.0 * vect[n - 1] + 17.0 * vect[n]);
//}

void printProgress(double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}
