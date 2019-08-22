//
// Created by Emlyn Graham on 9/08/19.
// Includes all custom global functions, structures and type definitions used in this library.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <complex>
#include <numeric>
#include <functional>
#include <algorithm>

#include "eigen/Eigen/Dense"


#ifndef EXTRAS_H
#define EXTRAS_H

using namespace Eigen;
using namespace std;

/// Typedefs

typedef ArrayXd dVec;
typedef ArrayXcd cVec;

typedef std::complex<double> cd;
typedef std::vector<int> intVec;
typedef std::vector<double> doubleVec;
typedef std::vector<std::complex<double>> complexVec;
typedef std::vector<complexVec> complexVecVec;
typedef std::vector<doubleVec> doubleVecVec;
typedef std::vector<std::pair<double, double> > doublePairVec;

/// Constants ///
extern std::complex<double> i;
extern double HBARC; // MeV fm
extern double AMU;   // MeV/c^2
extern double ESQ;    // Electron charge in MeV fm

/// Templates

template<typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


/// Structs

struct expVec { std::complex<double> operator()(std::complex<double> d) const { return std::exp(d); }};


/// Distribution functions

double sine(const double &x, const double &N, const double &xmin, const double &xmax);
double gaussian(const double &x, const double &X0, const double &Sigma);
double asymmGaussian(const double &x, const double &X0, const double &Sigma);


/// Vector functions

doubleVec vectorComplexToDouble(const complexVec &a);
complexVec vectorDoubleToComplex(const doubleVec &a);
complexVec vectorMultiply(const complexVec &a, const complexVec &b);
doubleVec vectorMultiply(const doubleVec &a, const doubleVec &b);
complexVec vectorExp(const complexVec &a);
complexVec vectorScale(const complexVec &a, const std::complex<double> &b);
complexVec vectorScale(const complexVec &a, const double &b);
doubleVec vectorScale(const doubleVec &a, const double &b);
complexVec vectorScale(const doubleVec &a, const std::complex<double> &b);
complexVec vectorAdd(const complexVec &a, const complexVec &b);
complexVec vectorAdd(const complexVec &a, const doubleVec &b);
complexVec vectorAdd(const doubleVec &a, const complexVec &b);
complexVec vectorSubtract(const complexVec &a, const complexVec &b);
complexVec vectorSubtract(const complexVec &a, const doubleVec &b);
complexVec vectorSubtract(const doubleVec &a, const complexVec &b);
doubleVec fourierComplexToDouble(const complexVec &cvector);
dVec fourierComplexToDouble(cVec &cvector);
complexVec fourierDoubleToComplex(const doubleVec &dvector);
cVec fourierDoubleToComplex(dVec &dvector);
double vectorSimpsonIntegrate(const doubleVec &vect, const double &h, const int &n);
double vectorSimpsonIntegrate(dVec &vect, const double &h, const int &n);
double vectorTrapezoidIntegrate(const doubleVec &vect, const double &h, const int &n);
double vectorTrapezoidIntegrate(dVec &vect, const double &h, const int &n);

#endif //EXTRAS_H
