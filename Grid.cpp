#include "Grid.h"
#include <iostream>
#include <vector>
#include <cmath>
typedef std::vector<double> double_vec;

Grid::Grid(int nstep, double xmin, double xmax, double kscale) {
  NStep = nstep;
  NPoint = NStep + 1;
  XMin = xmin;
  XMax = xmax;
  XStep = (XMax - XMin)/NStep;
  KScale = kscale;
  KStep = 2*M_PI/(NStep*XStep*KScale);
  KMin = -1*NStep*KStep/2.0;
  KMax = NStep*KStep/2.0;
  size_t size = NPoint;
  double_vec xarray(size);
  double_vec karray(size);
  double_vec earray(size);
  for(int i=0; i<NPoint; ++i){
    xarray[i] = i*XStep+XMin;
    karray[i] = i*KStep+KMin;
    earray[i] = pow((hbarc*karray[i]),2)/(2.0*amu);
  }
  X = xarray;
  K = karray;
  E = earray;
}

Grid::~Grid() {
  std::cout << "Object is being deleted" << std::endl;
}

void Grid::TestFcn() {
  std::cout << "Test Test" << std::endl;
}