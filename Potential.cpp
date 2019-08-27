#include <iostream>
#include <vector>

#include "Potential.h"

Potential::Potential(const Grid &object) : grid(1, 0.0, 1.0, 1.0) {
  grid = object;
  V.resize(grid.nPoint, Eigen::NoChange);
}

Potential::~Potential() {
  std::cout << "Potential deleted" << std::endl;
}

void Potential::test() {
  std::cout << "Test Potential" << std::endl;
  std::cout << "Potential is: " << V << std::endl;
}

/// Getters

dArray Potential::getReal() {
  return V.real();
}

dArray Potential::getImag() {
  return V.imag();
}

dArray Potential::getAbs() {
  return V.abs();
}

void Potential::initZero() {
  V.setZero(grid.nPoint);
}

/// Add to potential

void Potential::addConstant(const cd &c, const double &xmin, const double &xmax) {
  for (int j = 0; j < grid.nPoint; ++j) {
    if (grid.x(j) >= xmin and grid.x(j) <= xmax) {
      V(j) = V(j) + c;
    }
  }
}

void Potential::addParabolic(const double &xCentre, const cd &c) {
  V += c*(grid.x - xCentre).square();
}

void Potential::addQuartic(const double &xCentre, const cd &c) {
  V += c*(grid.x - xCentre).pow(4.0);
}

void Potential::addGaussian(const double &xCentre, const cd &height, const cd &sigma) {
  V += height*exp(-pow(grid.x - xCentre, 2.0)/(2.0*sigma*sigma));
}

void Potential::addWoodsSaxon(const double &xCentre,
                              const double &height,
                              const double &xSize,
                              const double &diffuseness) {
  V += 1+exp(((grid.x-xCentre).abs() - xSize)/diffuseness);
}

void Potential::addCoulomb(const double &Z1Z2, const double &xCentre, const double &xSize) {
  for (int j = 0; j < grid.nPoint; ++j) {
    if (std::abs(grid.x(j) - xCentre) < xSize) {
      V(j) = V(j) + cd(Z1Z2 * ESQ * (3.0 - pow(std::abs(grid.x(j) - xCentre) / xSize, 2.0) / (2.0 * xSize)), 0.0);
    } else {
      V(j) = V(j) + cd(((Z1Z2 * ESQ) / (std::abs(grid.x(j) - xCentre))), 0.0);
    }
  }
}

