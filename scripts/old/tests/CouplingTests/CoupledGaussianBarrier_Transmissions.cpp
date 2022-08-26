//
// Created by Emlyn Graham on 22/09/19.
// Tests error differences for various time and position steps
//
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <fstream>
#include <cmath>
#include "omp.h"
#include <string>

#include "Grid.h"
#include "Wavefunction.h"
#include "Potential.h"
#include "System.h"
#include "Plotter.h"

//#define EIGEN_RUNTIME_NO_MALLOC
using namespace Eigen;
using namespace std;

int main() {
  omp_set_num_threads(4);
  Eigen::setNbThreads(4);
//  Eigen::internal::set_is_malloc_allowed(false);
/// TODO: Fix System so you can add potentials and wavefunctions freely
/// Currently need to add wavefunctions first

  // First case in Dasso et al. 1983 'Channel coupling effects in heavy ion fusion reactions'

  string name = "DassoFig2";

  int time = 2900;
//  double timestep = 0.1;
  double timestep = 0.01;

  std::vector<double> Fs = {0.0, 2.0};

//  unsigned int sizeN = 2047;
  unsigned int sizeN = 16383;
  double xmin = -800.0;
  double xmax = 800.0;
  double kscale = 1.0;
  Grid grid(sizeN, xmin, xmax, kscale);
  for (int j = 0; j < Fs.size(); ++j) {
    std::vector<cd> initPsiK;
    double mu = 1.0;
    double V1 = 10.0;
    double V2 = 10.0;
    double sigma1 = 6.0;
    double sigma2 = 6.0;
    double F = Fs[j]; // coupling potential amplitude
    std::ofstream Out("DassoFig2"+tostring(int(F))+".csv");

    Wavefunction ground(grid, mu);

    ground.initGaussian(-50.0, 5.0);
    ground.boostEnergy(V1);
    Wavefunction excited(grid, mu);
    excited.initZero();

    Potential U1(grid);
    U1.initZero();
    U1.addGaussian(0.0, V1, 6.0);

    Potential U2(grid);
    U2.initZero();
    U2.addGaussian(0.0, V2, 6.0);

    Potential VC(grid);
    VC.initZero();
    VC.addGaussian(0.0, F, 6.0);

    System system(ground, U1);
    system.addWavefunction(excited);
    system.addPotential(VC, 0, 1);
    system.addPotential(VC, 1, 0);
    system.addPotential(U2, 1, 1);

    system.updateK();
    dArray groundPsiK_init = system.wavefunctions[0].psiK.abs();
    dArray excitedPsiK_init = system.wavefunctions[1].psiK.abs();

    // CC Evolution
    system.initCC(timestep);
    system.evolveCC(int(time/timestep));
    system.updateFromCC();
//
//    Plotter plot(system);
//    plot.animateCC(int(time/timestep), 100, false, false, true);

    dArray T = ((system.wavefunctions[0].psiK).abs2()+(system.wavefunctions[1].psiK).abs2())/((groundPsiK_init).abs2()+(excitedPsiK_init).abs2());

    Out << "E" << "," << "InitGround" << "," << "InitExcited" << "," << "FinalGround" << "," << "FinalExcited" << "," << "T\n";
    for (int j = ground.grid.nPoint/2; j < ground.grid.nPoint; ++j) {
      Out << system.wavefunctions[0].E(j) << "," << groundPsiK_init(j) << "," << excitedPsiK_init(j) << "," << abs(system.wavefunctions[0].psiK)(j) << "," << abs(system.wavefunctions[1].psiK)(j) << "," <<  T(j) << "\n";
    }
  }
}