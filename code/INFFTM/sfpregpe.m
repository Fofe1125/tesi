function [psi,a,b,Lmbda2,Lmbda] = sfpregpe(sf)
% Function file: [psi,a,b,Lmbda2,Lmbda] = sfpregpe(sf)
%
% Description: Compute preliminary quantities for a superfluids simulation.
%
% Input variables:
%   sf:     structure defining the function psi, where
%           sf.pdb:    are the physical domain boundaries
%           sf.psipdb: is the complex 3d-array of psi in the physical domain
%           sf.mirror: is the row vector of mirroring flags true/false
%           sf.t:      is simulation time (not used here)
%
% Output variables:
%   psi:    complex nd-array of the initial solution in the
%           computational domain
%   a:      row vector of left-hand points of the computational domain
%   b:      row vector of right-hand points of the computational domain
%   Lmbda2: eigenvalues of nd Laplace operator. They are defined as
%           Lmbda{1}.^2+Lmbda{2}.^2+... If psihat are the Fourier
%           coefficients of a function, Lmbda2.*psihat are the Fourier
%           coefficients of the Laplace operator applied to that function
%   Lmbda:  cell of coefficients for computing the derivatives of
%           Fourier series. If psihat are the Fourier coefficients of a
%           function, Lmbda{d}.*psihat are the Fourier coefficients of first
%           derivative in direction d
eta = size(sf.psipdb);
dd = sum(eta>1);
a = sf.pdb(1:2:2*dd-1);
b = sf.pdb(2:2:2*dd);
N = eta;
N(1:dd) = 2.^sf.mirror(1:dd).*eta(1:dd);
idx = cell(1,dd);
for d = 1:dd
  idx{d} = 1:eta(d);
end
oidx = idx;
psi = complex(zeros(N));
psi(idx{:}) = sf.psipdb;
for d = 1:dd
  if (sf.mirror(d))
    idx{d} = 1:N(d);
    psi(idx{:}) = cat(d,psi(oidx{:}),flip(psi(oidx{:}),d));
    oidx = idx;
    b(d) = b(d)+(b(d)-a(d));
  end
  lambda{d} = 2*pi*1i*(-N(d)/2:N(d)/2-1)/(b(d)-a(d));
end
[Lmbda{1:dd}] = ndgrid(lambda{1:dd});
Lmbda2 = sum(cat(dd+1,Lmbda{1:dd}).^2,dd+1); % Laplace eigenvalues
