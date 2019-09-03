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

//#include "eigen/Eigen/Dense"
//#include "eigen/unsupported/Eigen/CXX11/Tensor"
//#include "Eigen/Dense"
//#include "unsupported/Eigen/CXX11/Tensor"
#include "eigen/Eigen/Dense"
#include "eigen/unsupported/Eigen/CXX11/Tensor"

#ifndef EXTRAS_H
#define EXTRAS_H

using namespace Eigen;
using namespace std;

/// Typedefs

typedef std::complex<double> cd;
typedef std::vector<int> intVec;
typedef std::vector<double> doubleVec;
typedef std::vector<std::complex<double>> complexVec;
//typedef std::vector<complexVec> complexVecVec;
typedef std::vector<doubleVec> doubleVecVec;
typedef std::vector<std::pair<double, double> > doublePairVec;

typedef ArrayXd dArray;
typedef ArrayXcd cdArray;
typedef VectorXd dVector;
typedef VectorXcd cdVector;
typedef MatrixXcd cdMatrix;
typedef MatrixXd dMatrix;
typedef Eigen::Tensor<double, 3> dMatrixTensor;
typedef Eigen::Tensor<cd, 3> cdMatrixTensor;
typedef Eigen::Tensor<double, 2> dVectorTensor;
typedef Eigen::Tensor<cd, 2> cdVectorTensor;
typedef std::vector<cdVector> cdArrayVector;
typedef std::vector<cdMatrix> cdMatrixVector;
typedef std::vector<dVector> dArrayVector;
typedef std::vector<dMatrix> dMatrixVector;

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

template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using  VectorType = Eigen::Matrix<T,Eigen::Dynamic, 1>;


template<typename Scalar,int rank, typename sizeType>
auto Tensor_to_Matrix(const Eigen::Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols)
{
  return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows,cols);
}

template<typename Scalar, typename sizeType>
auto Tensor_to_Vector(const Eigen::Tensor<Scalar,1> &tensor,const sizeType rows)
{
  return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows);
}

template<typename Scalar, typename... Dims>
auto Matrix_to_Tensor(const MatrixType<Scalar> &matrix, Dims... dims)
{
  constexpr int rank = sizeof... (Dims);
  return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), {dims...});
}
template<typename Scalar, typename sizeType>
auto Vector_to_Tensor(const VectorType<Scalar> &matrix, const sizeType rows)
{
  return Eigen::TensorMap<Eigen::Tensor<const Scalar, 1>>(matrix.data(), rows);
}

template<typename T>
ostream& operator<< (ostream& out, const vector<T>& v) {
  out << "{";
  size_t last = v.size() - 1;
  for(size_t i = 0; i < v.size(); ++i) {
    out << v[i];
    if (i != last)
      out << ", ";
  }
  out << "}";
  return out;
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
doubleVec fourierComplexToDouble(const complexVec &cdVector);
dArray fourierComplexToDouble(cdArray &cdVector);
complexVec fourierDoubleToComplex(const doubleVec &dvector);
cdArray fourierDoubleToComplex(dArray &dvector);
//double vectorSimpsonIntegrate(const doubleVec &vect, const double &h, const int &n);
//double vectorSimpsonIntegrate(dArray &vect, const double &h, const int &n);
double vectorTrapezoidIntegrate(doubleVec &vect, const double &h, const int &n);
double vectorTrapezoidIntegrate(dArray &vect, const double &h, const int &n);

#endif //EXTRAS_H
