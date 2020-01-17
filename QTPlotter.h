//
// Created by Emlyn Graham on 13/08/19.
//

#ifndef PLOTTER_H
#define PLOTTER_H

#include "System.h"
#include "gnuplot-iostream/gnuplot-iostream.h"
#include "madplotlib/Madplotlib.h"

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
  double pauseTime = 1;
  //Functions
  explicit Plotter(System sys);
  ~Plotter();
  void test();
  /// TODO: Figure out best way to do this
  void plot();
  void plotPsi();
  void animate(int nSteps, double stepSize, int evolveOrder, int updateRate, bool logY);
  void animateCC(int nSteps, int updateRate, bool logY, bool k, bool cc);
  void animatePsi(int nSteps, double stepSize, int evolveOrder);
};

#endif //PLOTTER_H
