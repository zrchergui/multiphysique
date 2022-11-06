# multiphysique
Stage EV.

In this internship. I have used DuMuX, as simulation software. DuMuX is an open-source, C++ simulator from the University of Stuttgart designed to simulate flows in porous media. DuMuX is a toolbox based on DUNE (Distributed and Unified Numeric Environment), which was developed to solve a wide range of problems related to porous media. 

To solve the governing equations of the model, DuMuX uses monolithic  or sequential methods. The sequential algorithm is based on reformulating the equations of multi-phase flow into one equation for
pressure on the one hand,  and equations for phase and transport on the other hand. The most used sequential model is the fractional flow formulation for two-phase flows which is usually implemented applying an Implicit Pressure and Explicit Saturation algorithm (IMPES). 
 With the sequential structure, we can use different methods of discretization for different equations, the standard method used being a cell-centered finite-volume method. 
 
 The original module used in this internship is available at https://git.iws.uni-stuttgart.de/dumux-pub/Becker2018a  

 It contains two folders; one for each simulation type: energyStorage (full-dimensional model) and modelCoupling (adaptive model). Each folder contains several files, each having a specific functionality
main.cc: initializes parameters relative to the management of the grid geometry, the spatial parameters, the  boundary conditions, the solution vectors and the time parameters.
param.hh: creates tag types, calls other files, such as spatialparameters.hh  (cf. below) defines the discretization method used, calls the constitutive relations and the properties of the system and the state of the studied fluids, defines  properties such as the position  of the injectors, the injection rates and sets the initial conditions on the pressure.
spatialparam.hh: contains the parameters specific to the model studied such as porosity,  permeability, entry pressure, residual saturations, the $\gamma$-parameter used in the Brooks-Corey equations  for relative permeability and capillary pressure.
input: allows the user to enter the values of the execution parameters needed in the problem, such as inflow rate.
