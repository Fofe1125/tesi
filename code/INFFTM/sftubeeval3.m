function [Xi,rho] = sftubeeval3(sf,rhobar)
% Function file: [Xi,rho] = sftubeeval3(sf,rhobar)
%
% Description: Iteratively generate points such that rho <= rhobar(1:end).
%
% Input variables:
%   sf:     structure defining the function psi, where
%           sf.pdb:    are the physical domain boundaries
%           sf.psipdb: is the complex 3d-array of psi in the physical domain
%           sf.mirror: is the row vector of mirroring flags true/false
%           sf.t:      is the simulation time (not used here)
%   rhobar: decreasing sequence of density values delimiting vortex tubes
%
% Output variables:
%   Xi:     2d-array containing in row d the d-th components of the
%           final points corresponding to density rho <= rhobar(end)
%   rho:    row vector of density values corresponding to Xi
N = size(sf.psipdb);
rho = real(sf.psipdb).^2+imag(sf.psipdb).^2;
IND = find(rho <= rhobar(1));
lIND = length(IND);
[ix{1},ix{2},ix{3}] = ind2sub(N,IND);
[idx(:,1),idx(:,2),idx(:,3)] = ind2sub([3,3,3],(1:27)');
idx = idx-2;
% retrieve the grid in the physical domain
xcol = zeros(3,lIND*27);
x = cell(1,3);
h = zeros(1,3);
for d = 1:3
  x{d} = linspace(sf.pdb(2*d-1),sf.pdb(2*d),N(d)+1);
  x{d} = x{d}(1:N(d));
  h(d) = (sf.pdb(2*d)-sf.pdb(2*d-1))/N(d);
  temp = repmat(x{d}(ix{d}),27,1)+repmat(idx(:,d),1,lIND)*h(d)/3;
  xcol(d,:) = temp(:)';
end
[psihat,a,b] = sf2psihat(sf);
psixcol = infft3(psihat,a,b,xcol);
Xi = xcol;
rho = real(psixcol).^2+imag(psixcol).^2;
% refinements
kmax = length(rhobar);
for k = 2:kmax
  h = h/3;
  IND = find(rho <= rhobar(k));
  lIND = length(IND);
  Xi = Xi(:,IND);
  xcol = zeros(3,lIND*27);
  for d = 1:3
    temp = repmat(Xi(d,:),27,1)+repmat(idx(:,d),1,lIND)*h(d)/3;
    xcol(d,:) = temp(:)';
  end
  psixcol = infft3(psihat,a,b,xcol);
  Xi = xcol;
  rho = real(psixcol).^2+imag(psixcol).^2;
end
