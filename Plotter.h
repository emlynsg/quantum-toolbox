//
// Created by Emlyn Graham on 13/08/19.
//

#ifndef PLOTTER_H
#define PLOTTER_H

#include "System.h"
#include "gnuplot-iostream/gnuplot-iostream.h"

// http://stackoverflow.com/a/1658429
#ifdef _WIN32
#include <windows.h>
	inline void pause(unsigned millis) {
		::Sleep(millis);
	}
#else
#include <unistd.h>
inline void pause(unsigned millis) {
  ::usleep(millis * 1000);
}
#endif


class Plotter {
 public:
  // Attributes
  System system;
  int nColumns;
  bool showPsi;
  bool showPotential;
  bool showEnergy;
  bool showNorm;
  bool showAvgX;
  bool showPsiK;
  double pauseTime = 0.01;
  //Functions
  Plotter(System sys, const int &nColumns=2, const bool &showPsi=true,
          const bool &showPotential=false, const bool &showEnergy=false,
          const bool &showNorm=false, const bool &showAvgX=false, const bool &showPsiK=false);
  ~Plotter();
  void test();
  /// TODO: Figure out best way to do this
  void plot();
  void plotPsi();
  void animate(int nSteps, double stepSize, int evolveOrder);
  void animatePsi(int nSteps, double stepSize, int evolveOrder);
  void setPlotStyle(Gnuplot &g, int stylenum=0);
};

#endif //PLOTTER_H
