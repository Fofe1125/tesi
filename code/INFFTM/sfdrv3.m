% Description: Script to perform a 3d superfluids simulation.

%% Parameters defining spatial discretization in the physical domain
% physical domain boundaries
pdb = [-20,20,-20,20,-20,20];
% number (even) of grid points in pdb
eta = [40,40,40];
%% parameters defining the initial condition
% each row is a point in a vortex core
points = [2,0,0; -2,0,0];
% each row is the direction of a vortex core
vectors = [0,1,0; 0,0,1];
% if mirror(d), mirror in d direction
mirror = [true,true,true];
% anonymous function for the initial condition
ic = @(x) sfic3(x,pdb,mirror,points,vectors);
%% Parameters defining time discretization
% final time
T = 20;
% number of time steps
nsteps = 200;
%% Parameters defining output
% number of output files besides initial condition (it must divide nsteps)
nfiles = 10;
% root file name for output files
rootfilename = 'data';
%% end of user defined parameters
sf = sfrun(pdb,eta,ic,T,nsteps,nfiles,rootfilename);
