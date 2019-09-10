import matplotlib;
matplotlib.use("TkAgg");

import math;
import cmath;
import scipy.integrate;
import numpy;
import copy;
from matplotlib import pyplot as plt;
from matplotlib import animation;
import time;
import sys;

numpy.seterr(over="raise");

             # Times are in units of fm/c
hbarc=197.3; # MeV fm
amu=931.5;   # MeV/c^2
esq=1.44;    # Electron charge in MeV fm
             
# Grid class
class Grid :

    # This should define both spatial and momentum grids...g
    def __init__(self,NStep,XMin,XMax,KScale=1.0) :
        # Check number of steps on the grid
        self.NStep = NStep;
        self.NPoint = NStep+1;
        # Positions
        self.XMin = XMin;
        self.XMax = XMax;
        self.XStep = (self.XMax-self.XMin)/self.NStep;
        self.X = numpy.array([],dtype=float);
        # Momenta
        self.KScale=KScale;
        self.KStep = 2*math.pi/(self.NStep*self.XStep)*KScale;
        self.KMin = -self.NStep*self.KStep/2;
        self.KMax = self.NStep*self.KStep/2;
        self.K = numpy.array([],dtype=float);
	self.E = numpy.array([],dtype=float);
	# Set all the values
        for Index in range(0,self.NPoint,1) :
            self.X=numpy.append(self.X,Index*self.XStep + self.XMin);
            self.K=numpy.append(self.K,Index*self.KStep + self.KMin);
	    self.E=numpy.append(self.E,(hbarc*self.K[Index])**2/(2*amu));
        

# Wave function class
class WaveFunction :

    # Create wave function and initialize state to zero
    def __init__(self,Grid,ReducedMass=1) :
        global amu;
        self.Grid = Grid;
        self.Psi = numpy.zeros(self.Grid.NPoint,dtype="complex");
        self.PsiK = numpy.zeros(self.Grid.NPoint,dtype="complex");
        self.ReducedMass = ReducedMass * amu;

    def Normalise(self) :
        NormConstant = cmath.sqrt(Overlap(self.Grid,self.Psi,self.Psi));
        for Index in range(self.Grid.NPoint) :
            self.Psi[Index] = self.Psi[Index] / NormConstant;

    def CalculateNorm(self) :
        return Overlap(self.Grid,self.Psi,self.Psi);

    def CalculateNormInRegion(self,XMin,XMax) :
        PsiPart = numpy.zeros(self.Grid.NPoint,dtype="complex");
        for Index in range(self.Grid.NPoint) :
            if self.Grid.X[Index] > XMin and self.Grid.X[Index] < XMax :
                PsiPart[Index] = self.Psi[Index];
            else :
                PsiPart[Index] = 0.0;
        return Overlap(self.Grid,PsiPart,PsiPart);
        
    def ComputePsiK(self) :
        PsiFFT=numpy.multiply(self.Psi,numpy.exp(-1j*self.Grid.KMin*self.Grid.X));
        self.PsiK = scipy.fft(PsiFFT)*numpy.exp(-1j*self.Grid.K*self.Grid.XMin)*self.Grid.XStep/cmath.sqrt(2*math.pi);

    def ComputePsi(self) :
        PsiKFFT=numpy.multiply(self.PsiK,numpy.exp(1j*self.Grid.K*self.Grid.XMin)*cmath.sqrt(2*math.pi)/self.Grid.XStep);
        self.Psi = scipy.ifft(PsiKFFT)*numpy.exp(1j*self.Grid.KMin*self.Grid.X);
        
    def InitialZero(self) :
        for Index in range(0,self.Grid.NPoint,1) :
            self.Psi[Index] = 0.0;
            
    def InitialGaussian(self,X0,Sigma) :
        for Index in range(0,self.Grid.NPoint,1) :
            self.Psi[Index] = cmath.exp(-(self.Grid.X[Index]-X0)**2/(2*Sigma*Sigma));
        self.Normalise();

    def BoostWaveNumber(self,WaveNumber) :
        for Index in range(0,self.Grid.NPoint,1) :
            self.Psi[Index] = self.Psi[Index] * cmath.exp(1j*WaveNumber*self.Grid.X[Index]);
        self.Normalise();

    def Boost(self,Energy) :
        global hbarc;
        WaveNumber = numpy.sign(Energy) * cmath.sqrt(2*self.ReducedMass*abs(Energy))/hbarc;
        print(WaveNumber, Energy, self.ReducedMass,cmath.sqrt(2*self.ReducedMass*abs(Energy)));
        self.BoostWaveNumber(WaveNumber);
    
    def Real(self) :
        Result = [];
        for Index in range(0,self.Grid.NPoint,1) :
            Result.append(self.Psi[Index].real);
        return Result;
    
    def Imaginary(self) :
        Result = [];
        for Index in range(0,self.Grid.NPoint,1) :
            Result.append(self.Psi[Index].imag);
        return Result;

    def Absolute(self) :
        Result = [];
        for Index in range(0,self.Grid.NPoint,1) :
            Result.append(abs(self.Psi[Index]));
        return Result;

    def AbsoluteSquared(self) :
        Result = [];
        for Index in range(0,self.Grid.NPoint,1) :
            Result.append(abs(self.Psi[Index])**2);
        return Result;
    
    def RealK(self) :
        Result = [];
        for Index in range(0,self.Grid.NPoint,1) :
            Result.append(self.PsiK[Index].real);
        return Result;
    
    def ImaginaryK(self) :
        Result = [];
        for Index in range(0,self.Grid.NPoint,1) :
            Result.append(self.PsiK[Index].imag);
        return Result;

    def AbsoluteK(self) :
        Result = [];
        for Index in range(0,self.Grid.NPoint,1) :
            Result.append(abs(self.PsiK[Index]));
        return Result;

    def AverageX(self) :
        OverlapPsi = self.Grid.X*self.Psi;
        return Overlap(self.Grid,OverlapPsi,self.Psi);
    
class Potential :

    def __init__(self,Grid) :
        self.Grid = Grid;
        self.V = numpy.zeros(self.Grid.NPoint,dtype="complex");
    
    def Real(self) :
        Result = [];
        for Index in range(0,self.Grid.NPoint,1) :
            Result.append(self.V[Index].real);
        return Result;
    
    def Imaginary(self) :
        Result = [];
        for Index in range(0,self.Grid.NPoint,1) :
            Result.append(self.V[Index].imag);
        return Result;

    def Reset(self) :
        for Index in range(0,self.Grid.NPoint,1) :
            self.V[Index]=0.0;

    def SetConstant(self,Value,XMin,XMax) :
        for Index in range(self.Grid.NPoint) :
            if self.Grid.X[Index] >= XMin and self.Grid.X[Index] <= XMax :
                self.V[Index] = Value;
            
    def AddConstant(self,Value,XMin,XMax) :
        for Index in range(self.Grid.NPoint) :
            if self.Grid.X[Index] >= XMin and self.Grid.X[Index] <= XMax :
                self.V[Index] = self.V[Index] + Value;

    def AddParabolic(self,XCenter,Constant) :
        for Index in range(self.Grid.NPoint) :
            self.V[Index] = self.V[Index] + Constant*(self.Grid.X[Index]-XCenter)**2;

    def AddQuartic(self,XCenter,Constant) :
        for Index in range(self.Grid.NPoint) :
            self.V[Index] = self.V[Index] + Constant*(self.Grid.X[Index]-XCenter)**4;
            
    def AddGaussian(self,X0,Height,Sigma) :
        for Index in range(self.Grid.NPoint) :
            self.V[Index] = self.V[Index] + Height*cmath.exp(-(self.Grid.X[Index]-X0)**2/(2*Sigma*Sigma));

    def AddWoodsSaxon(self,Height,XCenter,XSize,Diffuseness) :
        for Index in range(self.Grid.NPoint) :
            if abs(self.Grid.X[Index]-XCenter) < 4*XSize :
                self.V[Index] = self.V[Index] + Height/(1+cmath.exp((abs(self.Grid.X[Index]-XCenter)-XSize)/Diffuseness));

    def AddCoulombSphere(self,Z1Z2,XCenter,XSize) :
        global esq;
        for Index in range(self.Grid.NPoint) :
            if abs(self.Grid.X[Index]-XCenter) < XSize :
                self.V[Index] = self.V[Index] + Z1Z2 * esq * (3.0 - (abs(self.Grid.X[Index]-XCenter)/XSize)**2 ) / (2*XSize);
            else :
                self.V[Index] = self.V[Index] + Z1Z2 * esq / abs(self.Grid.X[Index]-XCenter);
        
            
class System :

    def __init__(self,WF,Pot) :
        # Set the wave function and potential
        self.WF=WF;
        self.Pot=Pot;
        # Now set the arrays for energy tracking
        self.T=numpy.array([],dtype="complex");
        self.T=numpy.append(self.T,0.0);
        self.E=numpy.array([],dtype="float");
        self.E=numpy.append(self.E,self.Energy());
        self.Norm=numpy.array([],dtype="float");
        self.Norm=numpy.append(self.Norm,self.WF.CalculateNorm());
        self.AverageX=numpy.array([],dtype="float");
        self.AverageX=numpy.append(self.AverageX,self.WF.AverageX());

    # Expand the time evolution operator explicitly
    def Evolve(self,TimeStep,MaxOrder) :
        global hbarc;
        A = (hbarc/self.WF.Grid.XStep)**2 / (2*self.WF.ReducedMass);
        # Loop over all orders of the Taylor expansion of the evolution operator
        PsiPart = numpy.copy(self.WF.Psi);
        PsiTemp = numpy.copy(self.WF.Psi);
        # Now do the higher order terms
        for Order in range(1,MaxOrder+1,1) :
            # Do a fast numpy roll for all points
            PsiPart = numpy.add(A*(2*PsiTemp - numpy.roll(PsiTemp,-1) - numpy.roll(PsiTemp,1)),self.Pot.V*PsiTemp);
            # Now correct the edge terms
            PsiPart[0]=A*(2*PsiTemp[0] - PsiTemp[1]) + self.Pot.V[0]*PsiTemp[0];
            PsiPart[-1]=A*(2*PsiTemp[-1] - PsiTemp[-2]) + self.Pot.V[-1]*PsiTemp[-1];
            PsiTemp = -1*complex(0,1.0) * PsiPart * TimeStep / ( Order * hbarc );
            self.WF.Psi = numpy.add(self.WF.Psi,PsiTemp);
  
    def Log(self,TimeStep) :
        self.T=numpy.append(self.T,self.T[-1]+abs(TimeStep));
        self.E=numpy.append(self.E,self.Energy());
        self.Norm=numpy.append(self.Norm,self.WF.CalculateNorm());
        self.AverageX=numpy.append(self.AverageX,self.WF.AverageX());
        
    def Energy(self) :
        global hbarc;
        A = (hbarc/self.WF.Grid.XStep)**2 / (12*2*self.WF.ReducedMass);
        #OverlapPsi = numpy.add(A*(2*self.WF.Psi - numpy.roll(self.WF.Psi,-1) - numpy.roll(self.WF.Psi,1)),self.Pot.V*self.WF.Psi);
        OverlapPsi = numpy.add(A*(30*self.WF.Psi - 16*numpy.roll(self.WF.Psi,-1) - 16*numpy.roll(self.WF.Psi,1) + numpy.roll(self.WF.Psi,-2) + numpy.roll(self.WF.Psi,2)  ),self.Pot.V*self.WF.Psi);
        OverlapPsi[0] = 0.0;
        OverlapPsi[self.WF.Grid.NPoint-1] = 0.0;
        return Overlap(self.WF.Grid,OverlapPsi,self.WF.Psi);

class Plotter :

    def __init__(self,System,NumCol=2,ShowPot=False,ShowWF=True,ShowE=False,ShowNorm=False,ShowAverageX=False,ShowWFK=False) :
    
        # Create the plots prior to the evolution
        self.System=System;
        self.ShowPot=ShowPot;
        self.ShowWF=ShowWF;
        self.ShowE=ShowE;
        self.ShowNorm=ShowNorm;
        self.ShowAverageX=ShowAverageX;
        self.ShowWFK=ShowWFK;
        self.NumCol=NumCol;
        self.NPlots=0;
        self.NPointMax=2000;
        self.NPointAlways=2;
        
        # Create the figure
        self.Figure=plt.figure();
        self.AddLines();
        self.Figure.canvas.mpl_connect('resize_event',self.SaveBackgrounds);

        # Finally show the plot
        plt.show(block=False);
        
    def AddLines(self) :
        
        self.NumPlot=self.ShowPot+self.ShowWF+self.ShowE+self.ShowNorm+self.ShowAverageX+self.ShowWFK;
        self.NumRow=math.ceil(float(self.NumPlot)/self.NumCol);
        N=0;
        if self.NumPlot == 1 and self.ShowWF :
            self.NWFPlot=1;
        elif self.NumRow>0 :
            self.NWFPlot=self.NumCol+1;
        else :
            self.NWFPlot=2;

        # Potential
        if self.ShowPot :
            N=N+1;
            self.AxisPot=self.Figure.add_subplot(self.NumRow,self.NumCol,N);
            self.AxisPot.set_xlim([self.System.Pot.Grid.XMin,self.System.Pot.Grid.XMax]);
            self.LinePotReal,=self.AxisPot.plot([],[]);
            self.LinePotImag,=self.AxisPot.plot([],[]);
            self.AxisPot.set_xlabel("x [fm]");
            self.AxisPot.set_ylabel("Potential [MeV]");

        # Wave function at current time step
        if self.ShowWF :
            self.AxisWF=self.Figure.add_subplot(self.NumRow,self.NumCol,self.NWFPlot);
            self.AxisWF.set_xlim([self.System.WF.Grid.XMin,self.System.WF.Grid.XMax]);
            self.LineWFReal,=self.AxisWF.plot([],[]);
            self.LineWFImag,=self.AxisWF.plot([],[]);
            self.LineWFAbs, =self.AxisWF.plot([],[]);
            self.AxisWF.set_xlabel("x [fm]");
            self.AxisWF.set_ylabel("Psi(x)");

        # Wave function at current time step in momentum space
        if self.ShowWFK :
            N=N+1;
            if N==self.NWFPlot :
                N=N+1;
            self.AxisWFK=self.Figure.add_subplot(self.NumRow,self.NumCol,N);
            self.AxisWFK.set_xlim([self.System.WF.Grid.KMin,self.System.WF.Grid.KMax]);
            self.LineWFKReal,=self.AxisWFK.plot([],[]);
            self.LineWFKImag,=self.AxisWFK.plot([],[]);
            self.LineWFKAbs, =self.AxisWFK.plot([],[]);
            self.AxisWFK.set_xlabel("k [1/fm]");
            self.AxisWFK.set_ylabel("Psi(k)");

        # Energy as a function of time
        if self.ShowE :
            N=N+1;
            if N==self.NWFPlot :
                N=N+1;
            self.AxisE=self.Figure.add_subplot(self.NumRow,self.NumCol,N);
            self.AxisE.set_xlabel("Time [fm/c]");
            self.AxisE.set_ylabel("Energy [MeV]");
            self.LineE,=self.AxisE.plot(self.System.T,self.System.E);
            
        # Norm as a function of time
        if self.ShowNorm :
            N=N+1;
            if N==self.NWFPlot :
                N=N+1;
            self.AxisNorm=self.Figure.add_subplot(self.NumRow,self.NumCol,N);
            self.AxisNorm.set_xlabel("Time [fm/c]");
            self.AxisNorm.set_ylabel("Norm");
            self.LineNorm,=self.AxisNorm.plot(self.System.T,self.System.Norm);

        # Average position as a function of time
        if self.ShowAverageX :
            N=N+1;
            if N==self.NWFPlot :
                N=N+1;
            self.AxisAverageX=self.Figure.add_subplot(self.NumRow,self.NumCol,N);
            self.AxisAverageX.set_ylim([self.System.WF.Grid.X[0],self.System.WF.Grid.X[len(self.System.WF.Grid.X)-1]]);
            self.AxisAverageX.set_xlabel("Time [fm/c]");
            self.AxisAverageX.set_ylabel("Average X [fm]");
            self.LineAverageX,=self.AxisAverageX.plot(self.System.T,self.System.AverageX);

    def ClearLines(self) :
        if self.ShowPot :
            self.LinePotReal.set_xdata([]);
            self.LinePotReal.set_ydata([]);
            self.LinePotImag.set_xdata([]);
            self.LinePotImag.set_ydata([]);
        if self.ShowWF :
            self.LineWFReal.set_xdata([]);
            self.LineWFReal.set_ydata([]);
            self.LineWFImag.set_xdata([]);
            self.LineWFImag.set_ydata([]);
            self.LineWFAbs.set_xdata([]);
            self.LineWFAbs.set_ydata([]);
        if self.ShowWFK :
            self.LineWFKReal.set_xdata([]);
            self.LineWFKReal.set_ydata([]);
            self.LineWFKImag.set_xdata([]);
            self.LineWFKImag.set_ydata([]);
            self.LineWFKAbs.set_xdata([]);
            self.LineWFKAbs.set_ydata([]);
        if self.ShowE :
            self.LineE.set_xdata([]);
            self.LineE.set_ydata([]);
        if self.ShowNorm :
            self.LineNorm.set_xdata([]);
            self.LineNorm.set_ydata([]);
        if self.ShowAverageX :
            self.LineAverageX.set_xdata([]);
            self.LineAverageX.set_ydata([]);

    def Plot(self) :
        self.NPlots+=1;
        #if self.NPlots % 10 == 0 : self.SaveBackgrounds();
        if self.ShowPot :
            self.Figure.canvas.restore_region(self.BackgroundPot);
            self.LinePotReal.set_xdata(self.System.Pot.Grid.X);
            self.LinePotImag.set_xdata(self.System.Pot.Grid.X);
            self.LinePotReal.set_ydata(self.System.Pot.Real());
            self.LinePotImag.set_ydata(self.System.Pot.Imaginary());
            self.AxisPot.draw_artist(self.LinePotReal);
            self.AxisPot.draw_artist(self.LinePotImag);
            self.Figure.canvas.blit(self.AxisPot.bbox);
        if self.ShowWF :
            self.Figure.canvas.restore_region(self.BackgroundWF);
            self.LineWFReal.set_xdata(self.System.WF.Grid.X);
            self.LineWFImag.set_xdata(self.System.WF.Grid.X);
            self.LineWFAbs.set_xdata(self.System.WF.Grid.X);
            self.LineWFReal.set_ydata(self.System.WF.Real());
            self.LineWFImag.set_ydata(self.System.WF.Imaginary());
            self.LineWFAbs.set_ydata(self.System.WF.Absolute());
            self.AxisWF.draw_artist(self.LineWFReal);
            self.AxisWF.draw_artist(self.LineWFImag);
            self.AxisWF.draw_artist(self.LineWFAbs);
            self.Figure.canvas.blit(self.AxisWF.bbox);
        if self.ShowWFK :
            self.Figure.canvas.restore_region(self.BackgroundWFK);
            self.LineWFKReal.set_xdata(self.System.WF.Grid.K);
            self.LineWFKImag.set_xdata(self.System.WF.Grid.K);
            self.LineWFKAbs.set_xdata(self.System.WF.Grid.K);
            self.LineWFKReal.set_ydata(self.System.WF.RealK());
            self.LineWFKImag.set_ydata(self.System.WF.ImaginaryK());
            self.LineWFKAbs.set_ydata(self.System.WF.AbsoluteK());
            self.AxisWFK.draw_artist(self.LineWFKReal);
            self.AxisWFK.draw_artist(self.LineWFKImag);
            self.AxisWFK.draw_artist(self.LineWFKAbs);
            self.Figure.canvas.blit(self.AxisWFK.bbox);
        # Here we only plot a fraction of the recorded points (totalling
        # NPointMax) plus the last NPointAlways points.  This should keep the
        # data manageable for long calculations whilst still allowing the plot
        # to keep being updated.
        NSkip=int(math.floor(float(len(self.System.T))/self.NPointMax)+1);
        if self.ShowE :
            self.Figure.canvas.restore_region(self.BackgroundE);
            self.LineE.set_xdata(numpy.append(self.System.T[:-self.NPointAlways:NSkip],self.System.T[-self.NPointAlways:]));
            self.LineE.set_ydata(numpy.append(self.System.E[:-self.NPointAlways:NSkip],self.System.E[-self.NPointAlways:]));
            self.AxisE.draw_artist(self.LineE);
            self.Figure.canvas.blit(self.AxisE.bbox);
        if self.ShowNorm :
            self.Figure.canvas.restore_region(self.BackgroundNorm);
            self.LineNorm.set_xdata(numpy.append(self.System.T[:-self.NPointAlways:NSkip],self.System.T[-self.NPointAlways:]));
            self.LineNorm.set_ydata(numpy.append(self.System.Norm[:-self.NPointAlways:NSkip],self.System.Norm[-self.NPointAlways:]));
            self.AxisNorm.draw_artist(self.LineNorm);
            self.Figure.canvas.blit(self.AxisNorm.bbox);
        if self.ShowAverageX :
            self.Figure.canvas.restore_region(self.BackgroundAverageX);
            self.LineAverageX.set_xdata(numpy.append(self.System.T[:-self.NPointAlways:NSkip],self.System.T[-self.NPointAlways:]));
            self.LineAverageX.set_ydata(numpy.append(self.System.AverageX[:-self.NPointAlways:NSkip],self.System.AverageX[-self.NPointAlways:]));
            self.AxisAverageX.draw_artist(self.LineAverageX);
            self.Figure.canvas.blit(self.AxisAverageX.bbox);
        self.Figure.canvas.flush_events();
        
    def Autoscale(self) :
        if self.ShowPot :
            self.AxisPot.relim();
            self.AxisPot.autoscale();
        if self.ShowWF :
            self.AxisWF.relim();
            self.AxisWF.autoscale();
        if self.ShowWFK :
            self.AxisWFK.relim();
            self.AxisWFK.autoscale();
        if self.ShowE :
            self.AxisE.relim();
            self.AxisE.autoscale();
        if self.ShowNorm :
            self.AxisNorm.relim();
            self.AxisNorm.autoscale();
        if self.ShowAverageX :
            self.AxisAverageX.relim();
            self.AxisAverageX.autoscale();
        self.Figure.canvas.draw()
        self.Figure.canvas.flush_events();
        self.SaveBackgrounds();

    def SaveBackgrounds(self,event=None) :
        self.ClearLines();
        self.Figure.canvas.draw();
        if self.ShowPot :
            self.BackgroundPot=self.Figure.canvas.copy_from_bbox(self.AxisPot.bbox);
        if self.ShowWF: 
            self.BackgroundWF=self.Figure.canvas.copy_from_bbox(self.AxisWF.bbox);
        if self.ShowWFK: 
            self.BackgroundWFK=self.Figure.canvas.copy_from_bbox(self.AxisWFK.bbox);
        if self.ShowE :
            self.BackgroundE=self.Figure.canvas.copy_from_bbox(self.AxisE.bbox);
        if self.ShowNorm :
            self.BackgroundNorm=self.Figure.canvas.copy_from_bbox(self.AxisNorm.bbox); 
        if self.ShowAverageX :
            self.BackgroundAverageX=self.Figure.canvas.copy_from_bbox(self.AxisAverageX.bbox);
        self.Plot();

    def SavePlot(self,Filename) :
        self.Figure.savefig(Filename);
    
def Overlap(Grid,Psi1,Psi2) :
    if not numpy.isfinite(Psi1).any() or not numpy.isfinite(Psi2).any() :
        print(numpy.isfinite(Psi1), numpy.isfinite(Psi2));
        print("NON FINITE VALUES FOUND");
        sys.exit();
    Integrand=Psi1*numpy.conj(Psi2);
    return scipy.integrate.simps(Integrand,Grid.X);

def Transmission(Grid,PsiKi,PsiKf) :
    rhoKi=PsiKi*numpy.conj(PsiKi);
    rhoKf=PsiKf*numpy.conj(PsiKf);
    return rhoKf/rhoKi;
