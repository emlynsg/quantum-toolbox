//  omp_set_num_threads(4);
//  Eigen::setNbThreads(4);
//  Eigen::internal::set_is_malloc_allowed(false);
/// TODO: Fix System so you can add potentials and wavefunctions freely
/// Currently need to add wavefunctions first

// Linear scheme defined by
// N: Number of channels
// n \in {1,...,N}: Current channel
// V_{n->n+1}: n*VC Coupling strength
// Q_n: n*2 Relative energy

// Parabolic potential with constant strength couplings
// Linear and doubling coupling strengths
// Start with wavefunction in centre and boost

std::vector<int> configs = {3,4,5};
#pragma omp parallel for
for (int k = 0; k < configs.size(); ++k) {
int N = configs[k];
double mu = 1.0;
double V0 = 100.0;
double couplingCentre = 0.0;
double F = 4.0;
double baseEpsilon = 2.0; // Standard multiplier for epsilon
std::vector<double> ns;
std::vector<double> vcs;
std::vector<double> epsilons;
for (int j = 0; j < N; ++j) {
ns.push_back(j);
vcs.push_back(j*F);
//      vcs.push_back(pow(2.0,j)*F);
epsilons.push_back(j*baseEpsilon);
}

//    int time = 8500;
int time = 5000;
//    double timestep = 1.0;
//    double timestep = 0.1;
double timestep = 0.1;
//    unsigned int sizeN = 4095;
//    unsigned int sizeN = 8191;
unsigned int sizeN = 16383;
//    unsigned int sizeN = 32767;
//    unsigned int sizeN = 65535;
//    unsigned int sizeN = 131071;
double xmin = -1000.0;
double xmax = 1000.0;
double kscale = 1.0;

Grid grid(sizeN, xmin, xmax, kscale);

Wavefunction ground(grid, mu);
ground.initGaussian(0.0, 2.5);
ground.boostEnergy(V0);

Potential U1(grid);
U1.initZero();
U1.addParabolic(0.0, 0.001);

System system(ground, U1);
for (int j = 1; j < N; ++j) {
system.addZeroWavefunction(mu, epsilons[j]);
}
for (int j = 1; j < N; ++j) {
system.addParabolicPotential(0.0, 0.001, j, j);
}
// All coupled
for (int j = 0; j < N; ++j) {
for (int l = 0; l < N; ++l) {
if (j>l){
double randval = rand()/5.0;
system.addConstantPotential(randval, xmin, xmax, j, l);
system.addConstantPotential(randval, xmin, xmax, l, j);
cout << "Coupling " << j << " " << l << " is " << randval << endl;
}
}
}

// CC EvolutionC_
system.initCC(timestep);
std::vector<std::vector<double>> timeData;
std::vector<double> timeEnergies = {0.5*V0, 0.9*V0, 1.0*V0, 1.2*V0, 2.0*V0};
auto start = std::chrono::high_resolution_clock::now();
system.evolveCC(int(time/timestep), timeEnergies, timeData);
auto finish = std::chrono::high_resolution_clock::now();

for (int m = 0; m < timeEnergies.size(); ++m) {
//      std::ofstream Out("N_"+tostring(N)+"_E_"+tostring(timeEnergies[m]/V0)+".csv");
std::ofstream Out("N_"+tostring(N)+"_E_"+tostring(timeEnergies[m])+".csv");
int len = timeData[0].size();
Out << "t" << "," << "T+R";
for (int l = 0; l < system.nChannel; ++l) {
Out << "," << "N_"+tostring(l);
}
Out << "\n";
for (int j = 0; j < len; ++j) {
Out << timeData[0][j] << "," << timeData[1+m*(system.nChannel+1)][j];
for (int l = 0; l < system.nChannel; ++l) {
Out << "," << timeData[1+m*(system.nChannel+1)+1+l][j];
}
Out << "\n";
}
}
system.updateFromCC();
//
//    Plotter plot(system);
//    plot.animateCC(int(time/timestep), 100, false, false, true);
std::ofstream Out2("N_"+tostring(N)+"ChannelsOverTime.csv");
Out2 << "t";
for (int l = 0; l < system.nChannel; ++l) {
Out2 << ",Norm_"+tostring(l);
Out2 << ",AvgX_"+tostring(l);
}
Out2 << "\n";
for (int l = 0; l < system.times[0].size(); ++l) {
Out2 << system.times[0][l];
for (int m = 0; m < system.nChannel; ++m) {
Out2 << "," << system.norms[m][l];
Out2 << "," << system.averages[m][l];
}
Out2 << "\n";
}

std::chrono::duration<double> elapsed = finish - start;
std::cout << "Elapsed time for "+tostring(N)+" couplings: " << elapsed.count() << " s\n";

//    dArray T = system.getTransmission();

//    std::vector<dArray> Rs;
//    for (int l = 0; l < N; ++l) {
//      Rs.push_back(system.getReflection(l));
//    }
//
//    std::ofstream Out("N_"+tostring(N)+".csv");
//    Out << "E,T,";
//    for (int l = 0; l < N; ++l) {
//      Out << "R"+tostring(l)+",";
//    }
//    Out << "R,";
//    for (int l = 0; l < N; ++l) {
//      Out << "N_R"+tostring(l)+",";
//    }
//    Out << "E/V0\n";
//    for (int j = 1+ground.grid.nPoint/2; j < ground.grid.nPoint; ++j) {
//      Out << system.wavefunctions[0].E(j) << "," << T(j) << ",";
//      double tot = 0.0;
//      for (int l = 0; l < N; ++l) {
//        double Rl = Rs[l](j);
//        Out << Rl << ",";
//        tot += Rl;
//      }
//      Out << tot << ",";
//      for (int l = 0; l < N; ++l) {
//        double Rl = Rs[l](j);
//        Out << Rl/(tot+T(j)) << ",";
//      }
//      Out << system.wavefunctions[0].E(j)/V0 << "\n";
//    }
}