//
// Created by emlyn on 25/02/2019.
//
#include "Grid.h"
#include <iostream>
#include <math.h>

using namespace std;
class Grid
{
public:
    // Attributes
    // (Issue: want to be able to create some of these as arrays given other inputs, need to figure this out)
    double NStep;
    double NPoint;
    double XMin;
    double XMax;
    double KScale;
    double KStep;
    double[] X;
    double[] K;
    double[] E;

    //Parametrized Constructor
    Grid(n_step,x_min,x_max,k_scale)
    {
        cout << "Parametrized Constructor called" << endl;
        NStep = n_step;
        XMin = x_min;
        XMax = x_max;
        KScale = k_scale;
        //Check number of steps on the grid
        NPoint = NStep+1.0;
        //Positions
        XStep = (XMax-XMin)/NStep;
        X = ;//Something - figure this out (numpy.array([],dtype=float))
        //Momenta
        KStep = 2.0*3.141592653589793238463/(NStep*XStep)*KScale;
        KMin = -NStep*KStep/2.0;
        KMax = NStep*KStep/2.0;
        K = ;//Something - figure this out (numpy.array([],dtype=float))
        E = ;//Something - figure this out (numpy.array([],dtype=float))
        //Set all the values
        for(i = 0; i < NPoint; i++)
        {
            X[i] = (X,i*XStep+XMin);
            K[i] = ;
            E[i] = ;
        }
        for Index in range(0,self.NPoint,1) :
        self.X=numpy.append(self.X,Index*self.XStep + self.XMin);
        self.K=numpy.append(self.K,Index*self.KStep + self.KMin);
        self.E=numpy.append(self.E,(hbarc*self.K[Index])**2/(2*amu));

    }
};