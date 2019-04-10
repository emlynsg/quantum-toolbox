#include <iostream>
#include <vector>
#include "Grid.h"
typedef std::vector<double> double_vec;

int main() {

  /// We need to first make some variables for our Grid class, and then put them into the constructor ///

  int sizeN = 10;
  double xmin = 0.0;
  double xmax = 1.0;
  double kscale = 1.0;


  size_t size = 10;
  double_vec gridarray(size);
  for(int i=0; i<size; ++i){
    gridarray[i] = i/10.0;
  }

  Grid gridObject(sizeN, xmin, xmax, kscale);

  /// Checking that the class object was instantiated properly ///

  gridObject.TestFcn();
  std::cout << gridObject.X[4] << std::endl;



  return 0;

}