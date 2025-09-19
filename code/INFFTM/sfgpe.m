function psi = sfgpe(psi,T,nsteps,Lmbda2)
% Function file: psi = sfgpe(psi,T,nsteps,Lmbda2)
%
% Description: Solve superfluid GPE by time splitting Fourier spectral method.
%
% Input variables:
%   psi:    complex nd-array of the initial solution in the
%           computational domain
%   T:      final computational time
%   nsteps: number of time steps
%   Lmbda2: eigenvalues of nd Laplace operator. They are defined as
%           Lmbda{1}.^2+Lmbda{2}.^2+... If psihat are the Fourier
%           coefficients of a function, Lmbda2.*psihat are the Fourier
%           coefficients of the Laplace operator applied to that function
%
% Output variables:
%   psi:    complex nd-array of solution at time T in the computational domain
k = T / nsteps;
expLambda2 = exp(k* 1i/2 * Lmbda2);
V = ones(size(psi)); % scalar potential
psi = exp(k/2 * (1i/2 * (V-real(psi).^2-imag(psi).^2))) .* psi;
for j = 1:nsteps-1
  psi = fftshift(fftn(psi));%*prod(sqrt(b-a)./N);
  psi = expLambda2 .* psi;
  psi = ifftn(ifftshift(psi));%/prod(sqrt(b-a)./N);
  psi = exp(k * (1i/2 * (V-real(psi).^2-imag(psi).^2))) .* psi;
end
psi = fftshift(fftn(psi));%*prod(sqrt(b-a)./N);
psi = expLambda2 .* psi;
psi = ifftn(ifftshift(psi));%/prod(sqrt(b-a)./N);
psi = exp(k/2 * (1i/2 * (V-real(psi).^2-imag(psi).^2))) .* psi;
