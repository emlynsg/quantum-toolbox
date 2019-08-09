#include "Grid.h"

Grid::Grid(int nstep, double xmin, double xmax, double kscale) {
  nStep = nstep;
  nPoint = nStep + 1;
  xMin = xmin;
  xMax = xmax;
  xStep = (xMax - xMin) / nStep;
  kScale = kscale;
  kStep = 2 * M_PI / (nStep * xStep * kScale);
  kMin = -1 * nStep * kStep / 2.0;
  kMax = nStep * kStep / 2.0;
  x.reserve(nPoint);
  k.reserve(nPoint);
  E.reserve(nPoint);
  for (int j = 0; j < nPoint; ++j) {
    x.push_back(j * xStep + xMin);
    k.push_back(j * kStep + kMin);
    E.push_back(pow((HBARC * k[j]), 2.0) / (2.0 * AMU));
  }
}

Grid::~Grid() {
  doubleVec().swap(x);
  doubleVec().swap(k);
  doubleVec().swap(E);
  std::cout << "Grid deleted" << std::endl;
}

void Grid::test() {
  std::cout << "Test Grid" << std::endl;
  std::cout << "Number of grid points: " << nPoint << std::endl;
}