//
// Created by Emlyn Graham on 13/08/19.
//

#ifndef PLOTTER_H
#define PLOTTER_H

#include "System.h"
#include "gnuplot-cpp/include/gnuplot.h"

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
  //Functions
  Plotter(System sys, const int &nColumns=2, const bool &showPsi=true,
          const bool &showPotential=false, const bool &showEnergy=false,
          const bool &showNorm=false, const bool &showAvgX=false, const bool &showPsiK=false);
  ~Plotter();
  void test();
  void plot();
  void setPlotStyle(GnuplotPipe &g, int stylenum=0);

};

#endif //PLOTTER_H
