#include "Wavefunction.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "Grid.h"
typedef std::vector<double> double_vec;
typedef std::vector< std::complex<double> > complex_vec;
typedef std::complex<double> complex;

Wavefunction::Wavefunction(const Grid& object, double ReducedMass) : grid(1,0.0,1.0,1.0) {
  grid = object;
  reduced_mass = ReducedMass*amu;
  for(int i=0; i<grid.n_point; ++i){
    psi.push_back(complex(0.0, 0.0));
    psi_k.push_back(complex(0.0, 0.0));
  }

}

Wavefunction::~Wavefunction() {
  std::cout << "Object is being deleted" << std::endl;
}

void Wavefunction::TestFcn() {
  std::cout << "Test Test" << std::endl;
}