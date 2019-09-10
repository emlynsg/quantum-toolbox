import qmtoybox as qm;

# Illustrates tunnelling through a square barrier

# Parameters of the time evolution
TimeStep=0.1
NTimeStep=100000;
EvolveOrder=20;

# A parameter to control how often the plots are updated
UpdateInterval=100;

# Establish a grid for the calculations
Grid = qm.Grid(4095,-150.0,150.0);

# Create a wave function and assign some values
WaveFunction = qm.WaveFunction(Grid);
WaveFunction.InitialGaussian(-100.0,10.0);
WaveFunction.Boost(80.0);

# Create a potential
Potential = qm.Potential(Grid);
Potential.AddConstant(75.0, -2.0, -1.0)

# Create the system
System = qm.System(WaveFunction,Potential);

# Create the plotter and set properties
Plotter = qm.Plotter(System,2,ShowPot=True,ShowWF=True,ShowWFK=True,ShowNorm=True,ShowE=True,ShowAverageX=True);
Plotter.AxisPot.set_ylim([-50.0,120.0]);
Plotter.AxisWF.set_ylim([-0.5,0.5]);
Plotter.AxisWFK.set_xlim([-5.0,5.0]);
Plotter.AxisWFK.set_ylim([-3.0,3.0]);
Plotter.AxisE.set_xlim([0,abs(TimeStep*NTimeStep)]);
Plotter.AxisE.set_ylim([0.1,1000]);
Plotter.AxisE.set_yscale('log');
Plotter.AxisNorm.set_xlim([0,abs(TimeStep*NTimeStep)]);
Plotter.AxisNorm.set_ylim([0.001,1.1]);
Plotter.AxisNorm.set_yscale('log');
Plotter.AxisAverageX.set_xlim([0,abs(TimeStep*NTimeStep)]);
Plotter.SaveBackgrounds();

# Loop in time to do the time evolution
for Index in xrange(NTimeStep):
    # Call System.Evolve to do a single step in the time evolution
    System.Evolve(TimeStep,EvolveOrder);
    # Update the plots every UpdateInterval time steps
    if Index % UpdateInterval == 0 :
        # Compute the Fourier transform of the state
        WaveFunction.ComputePsiK();
        # Log the energy, average position etc.
        System.Log(UpdateInterval*TimeStep);
        # Update the actual plot
        Plotter.Plot();

# Write the final state to a file
Plotter.SavePlot("squarebarrier.png");
