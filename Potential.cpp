//# define NDEBUG

#include <iostream>
#include <vector>

#include "Potential.h"

Potential::Potential(const Grid &object) : grid(1, 0.0, 1.0, 1.0) {
  grid = object;
  V.resize(grid.nPoint, 0.0);
}

Potential::~Potential() {
  complexVec().swap(V);
  std::cout << "Potential deleted" << std::endl;
}

void Potential::test() {
  std::cout << "Test Potential" << std::endl;
  std::cout << "First element: " << V[0] << std::endl;
}

/// Getters

doubleVec Potential::getReal() {
  doubleVec returnValue(grid.nPoint);
  std::transform(V.begin(), V.end(), returnValue.begin(), [](auto &elt) { return elt.real(); });
  return returnValue;
}

doubleVec Potential::getImag() {
  doubleVec returnValue(grid.nPoint);
  std::transform(V.begin(), V.end(), returnValue.begin(), [](auto &elt) { return elt.imag(); });
  return returnValue;
}

doubleVec Potential::getAbs() {
  doubleVec returnValue(grid.nPoint);
  std::transform(V.begin(), V.end(), returnValue.begin(), [](auto &elt) { return std::abs(elt); });
  return returnValue;
}

void Potential::initZero() {
  std::fill(V.begin(), V.end(), 0.0);
}

void Potential::initConstantInRegion(const double &c, const double &xmin, const double &xmax) {
  for (int j = 0; j < grid.nPoint; ++j) {
    if (grid.x[j] >= xmin and grid.x[j] <= xmax) {
      V[j] = complex(c, 0.0);
    }
  }
}

/// Add to potential

void Potential::addConstant(const double &c, const double &xmin, const double &xmax) {
  for (int j = 0; j < grid.nPoint; ++j) {
    if (grid.x[j] >= xmin and grid.x[j] <= xmax) {
      V[j] = V[j] + complex(c, 0.0);
    }
  }
}

void Potential::addParabolic(const double &xCenter, const double &c) {
  std::transform(V.begin(),
                 V.end(),
                 grid.x.begin(),
                 V.begin(),
                 [c, xCenter](auto &v, auto &x) { return v + complex((c * pow(x - xCenter, 2.0)), 0.0); });
}

void Potential::addQuartic(const double &xCenter, const double &c) {
  std::transform(V.begin(),
                 V.end(),
                 grid.x.begin(),
                 V.begin(),
                 [c, xCenter](auto &v, auto &x) { return v + complex((c * pow(x - xCenter, 4.0)), 0.0); });
}

void Potential::addGaussian(const double &xCenter, const double &height, const double &sigma) {
  std::transform(V.begin(),
                 V.end(),
                 grid.x.begin(),
                 V.begin(),
                 [xCenter, height, sigma](auto &v, auto &x) {
                   return v + complex(exp(-pow(x - xCenter, 2.0) / (2.0 * sigma * sigma)), 0.0);
                 });
}

void Potential::addWoodsSaxon(const double &xCenter,
                              const double &height,
                              const double &xSize,
                              const double &diffuseness) {
  std::transform(V.begin(),
                 V.end(),
                 grid.x.begin(),
                 V.begin(),
                 [xCenter, height, xSize, diffuseness](auto &v, auto &x) {
                   return v + complex((1 + exp((abs(x - xCenter) - xSize) / diffuseness)), 0.0);
                 });
}

void Potential::addCoulombSphere(const double &Z1Z2, const double &xCenter, const double &xSize) {
  for (int j = 0; j < grid.nPoint; ++j) {
    if (grid.x[j] - xCenter < xSize) {
      V[j] = V[j] + complex(Z1Z2 * ESQ * (3.0 - pow(std::abs(grid.x[j] - xCenter) / xSize, 2.0) / (2.0 * xSize)), 0.0);
    } else {
      V[j] = V[j] + complex(Z1Z2 * ESQ / (std::abs(grid.x[j] - xCenter)), 0.0);
    }
  }
}

