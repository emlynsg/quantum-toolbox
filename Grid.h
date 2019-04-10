#ifndef PRACTICECODE_GRID_H
#define PRACTICECODE_GRID_H
#include <vector>
typedef std::vector<double> double_vec;


class Grid {

 public:
  /// Constants ///
  double hbarc=197.3; /// MeV fm
  double amu=931.5;   /// MeV/c^2
  double esq=1.44;    /// Electron charge in MeV fm
  /// Number of steps on the grid ///
  int NStep;
  int NPoint;
  /// Grid positions ///
  double XMin;
  double XMax;
  double XStep;
  double_vec X;
  /// Momenta ///
  double KScale;
  double KStep;
  double KMin;
  double KMax;
  double_vec K;
  /// Energies ///
  double_vec E;
  Grid(int nstep, double xmin, double xmax, double kscale);
  ~Grid();
  void TestFcn();

};

#endif //PRACTICECODE_GRID_H
