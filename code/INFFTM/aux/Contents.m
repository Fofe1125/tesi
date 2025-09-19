% INFFTM: a MATLAB code for the evaluation of 3d Fourier series at
% arbitrary points and for the computation of the solution of the
% Gross-Pitaevskii equation for superfluids as a 3d Fourier series
%
% Auxiliary functions.
%
% fftshift   - Fix for bug #45207 in GNU Octave 4.0.0.
% flip       - flip.m function for GNU Octave < 4.0.0.
% ifftshift  - Fix for bug #45207 in GNU Octave 4.0.0.
% igridft2   - Evaluate 2-dimensional Fourier coeffs to a rectilinear grid.
% infft2     - Evaluate 2d-array of Fourier coefficients by NFFT.
% infft      - Evaluate 1d-array of Fourier coefficients by NFFT.
% mlloadnfft - Script to fix bug. no.961694 in certain Matlab versions.
% nfftpath   - Script for selecting NFFT usage.
%
% Authors:
% Marco Caliari and Simone Zuccher
