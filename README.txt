%===============================%
% Manhattan Lattice             %
% Addison Hartman               %
% Weld Lab at UCSB              %
% Last edited in October 2020   %
%===============================%

%===============================%
   Manhattan.m
%===============================%

The Matlab code file Manhattan.m allows you to run simulations of a BEC evolving through the Manhattan lattice. The two sections which are meant to be edited routinely for each simulation are:

- "File outputs": determine which plots and gifs are produced and the settings for the gif

- "Simulation parameters": change the size, mesh spacing, number of saved wave functions for the gif, ramp time, end time, and time step size

The rest of the sections should NOT be edited for a typical simulation. However, they need to be edited to improve the simulations. The next suggested edits are:

- make gif creation into a function for easier readability
- include interactions by adjusting the position propagator and changing the initial state to a Thomas-Fermi distribution
- make a momentum gif for comparison to in-lab measurements

%===============================%
   gifs folder
%===============================%

The gifs folder contains a multitude of gifs created with this code. They are numbered first by date in the form YYMMDD-HHMM, and then have a string of parameters for easy searching. Note that all files preceding 201002-0912 do not have the correct timing indicated.