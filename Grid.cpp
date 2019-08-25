#include "Grid.h"

using namespace Eigen;
using namespace std;

Grid::Grid(unsigned int nstep, double xmin, double xmax, double kscale) {
  nStep = nstep;
  nPoint = nStep + 1;
  xMin = xmin;
  xMax = xmax;
  xStep = (xMax - xMin) / nStep;
  kScale = kscale;
  kStep = 2.0 * M_PI / (nStep * xStep * kScale);
  kMin = -1.0 * nStep * kStep / 2.0;
  kMax = nStep * kStep / 2.0;
  x.resize(nPoint, Eigen::NoChange);
  k.resize(nPoint, Eigen::NoChange);
  E.resize(nPoint, Eigen::NoChange);
  x = ArrayXd::LinSpaced(nPoint, xMin, xMax);
  k = ArrayXd::LinSpaced(nPoint, kMin, kMax);
  E = ((HBARC*k).square())/(2.0*AMU);
}

Grid::~Grid() {
//  std::cout << "Grid deleted" << std::endl;
}

void Grid::test() {
  std::cout << "Test Grid" << std::endl;
  std::cout << "Grid is: " << x << std::endl;
}