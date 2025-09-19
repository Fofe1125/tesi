% INFFTM: a MATLAB code for the evaluation of 3d Fourier series at
% arbitrary points and for the computation of the solution of the
% Gross-Pitaevskii equation for superfluids as a 3d Fourier series
%
% Drivers and main functions.
%
% nfftpath    - Script for selecting NFFT usage.
% sfdrv3      - Script to perform a 3d superfluids simulation.
% evaldrv3    - Script to perform fast evaluation at grids or artibrary points.
% igridftn    - Evaluate nd-array of Fourier coeffs to a rectilinear grid.
% infft3      - Evaluate 3d-array of Fourier coefficients by NFFT.
% mkgrid      - Generate a rectilinear grid in the physical domain.
% ndcovlt     - Compute an n-dimensional linear transform of an nd tensor.
% plotiso3    - Plot an isosurface of the real 3d input data.
% sf2psihat   - Extract nd-array of Fourier coefficients from a sf structure.
% sf4pade     - Compute Pade' [8/8] approximation of a 2d vortex density.
% sfEm        - Compute mass and energy of a superfluid in the physical domain.
% sfgpe       - Solve superfluid GPE by time splitting Fourier spectral method.
% sfic3       - Compute the structure sf for superimposed 3d straight vortices.
% sfpregpe    - Compute preliminary quantities for a superfluids simulation.
% sfrun       - Run a complete superfluid simulation.
% sfsvl3      - Compute psi originating from a single straight 3d vortex.
% sftubeeval3 - Iteratively generate points such that rho <= rhobar(1:end).
% sfview3     - Show an isosurface of the density associated to a structure sf.
%
% aux         - Folder containing auxiliary functions.
%
% Authors:
% Marco Caliari and Simone Zuccher
