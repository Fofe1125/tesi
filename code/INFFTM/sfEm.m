function [mass,enrgy,Lmbda2] = sfEm(sf,Lmbda2)
% Function file: [mass,enrgy] = sfEm(sf)
% Function file: [mass,enrgy] = sfEm(sf,Lmbda2)
% Function file: [mass,enrgy,Lmbda2] = sfEm(sf)
%
% Description: Compute mass and energy of a superfluid in the physical domain.
%
% Input variables
%   sf:     structure defining the function psi, where
%           sf.pdb:    are the physical domain boundaries
%           sf.psipdb: is the complex 3d-array of psi in the physical domain
%           sf.mirror: is the row vector of mirroring flags true/false
%           sf.t:      is simulation time (not used here)
%   Lmbda2: (optional) eigenvalues of nd Laplace operator. They are defined as
%           Lmbda{1}.^2+Lmbda{2}.^2+... If psihat are the Fourier
%           coefficients of a function, Lmbda2.*psihat are the Fourier
%           coefficients of the Laplace operator applied to that function
%
% Output variables:
%   mass:   mass in the physical domain
%   enrgy:  energy in the physical domain
%   Lmbda2: (optional) eigenvalues of nd Laplace operator. They are defined as
%           Lmbda{1}.^2+Lmbda{2}.^2+... If psihat are the Fourier
%           coefficients of a function, Lmbda2.*psihat are the Fourier
%           coefficients of the Laplace operator applied to that function
[psihat,a,b,psi] = sf2psihat(sf);
N = size(psi);
dd = sum(N>1);
V = ones(prod(N),1); % potential
h = (b - a) ./ N;
E2 = prod(h)*sum((V-(real(psi(:)).^2+imag(psi(:)).^2)).^2)/4;
temp = real(psihat).^2+imag(psihat).^2;
mass = sum(temp(:))/2^sum(sf.mirror);
if (nargin == 1)
  for d = 1:dd
    lambda{d} = 2*pi*1i*(-N(d)/2:N(d)/2-1)/(b(d)-a(d));
  end
  [Lmbda{1:dd}] = ndgrid(lambda{1:dd});
  Lmbda2 = sum(cat(dd+1,Lmbda{1:dd}).^2,dd+1); % Laplace eigenvalues
end
E1 = -sum(Lmbda2(:).*temp(:))/2;
enrgy = (E1+E2)/2^sum(sf.mirror);
