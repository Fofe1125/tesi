function [psihat,a,b,psi] = sf2psihat(sf)
% Function file: [psihat,a,b,psi] = sf2psihat(sf)
%
% Description: Extract nd-array of Fourier coefficients from a sf structure.
%
% Input variables:
%   sf:     structure defining the function psi, where
%           sf.pdb:    are the physical domain boundaries
%           sf.psipdb: is the complex nd-array of psi in the physical domain
%           sf.mirror: is the row vector of mirroring flags true/false
%           sf.t:      is the simulation time (not used here)
%
% Output variables:
%   psihat: complex nd-array of Fourier coefficients of psi
%   a:      row vector of left-hand points of the computational domain
%   b:      row vector of right-hand points of the computational domain
%   psi:    complex nd-array of psi in the computational domain
eta = size(sf.psipdb);
dd = sum(eta>1);
a = sf.pdb(1:2:2*dd-1);
b = sf.pdb(2:2:2*dd);
N = eta;
N(1:dd) = 2.^sf.mirror(1:dd).*eta(1:dd);
psi = complex(zeros(N));
idx = cell(1,dd);
for d = 1:dd
  idx{d} = 1:eta(d);
end
oidx = idx;
psi(idx{:}) = sf.psipdb;
for d = 1:dd
  if (sf.mirror(d))
    idx{d} = 1:N(d);
    psi(idx{:}) = cat(d,psi(oidx{:}),flip(psi(oidx{:}),d));
    oidx = idx;
    b(d) = b(d)+(b(d)-a(d));
  end
end
psihat = fftshift(fftn(psi)) * prod(sqrt(b-a) ./ N(1:dd));
%!test
%! a = [-10,-12,-14];
%! b = [10,12,14];
%! eta = [10,12,14];
%! sf.mirror = true*ones(1,3);
%! N = eta;
%! N = 2.^sf.mirror.*eta;
%! psi = rand(eta);
%! sf.psipdb = psi;
%! sf.pdb = reshape([a;b],1,6);
%! for d = 1:3
%!   if (sf.mirror(d))
%!     psi = cat(d,psi,flip(psi,d));
%!     b(d) = b(d)+(b(d)-a(d));
%!   end
%! end
%! psihatref = fftshift(fftn(psi)) * prod(sqrt(b-a) ./ N);
%! psihat = sf2psihat(sf);
%! assert(psihatref(:),psihat(:))
%!test
%! a = [-10,-12];
%! b = [10,12];
%! eta = [10,12];
%! sf.mirror = true*ones(1,2);
%! N = eta;
%! N = 2.^sf.mirror.*eta;
%! psi = rand(eta);
%! sf.psipdb = psi;
%! sf.pdb = reshape([a;b],1,4);
%! for d = 1:2
%!   if (sf.mirror(d))
%!     psi = cat(d,psi,flip(psi,d));
%!     b(d) = b(d)+(b(d)-a(d));
%!   end
%! end
%! psihatref = fftshift(fft2(psi)) * prod(sqrt(b-a) ./ N);
%! psihat = sf2psihat(sf);
%! assert(psihatref(:),psihat(:))
%!test
%! a = -10;
%! b = 10;
%! eta = 10;
%! sf.mirror = true;
%! N = eta;
%! N = 2.^sf.mirror.*eta;
%! psi = rand(eta,1);
%! sf.psipdb = psi;
%! sf.pdb = reshape([a;b],1,2);
%! if (sf.mirror(1))
%!   psi = cat(1,psi,flip(psi,1));
%!   b(1) = b(1)+(b(1)-a(1));
%! end
%! psihatref = fftshift(fft(psi)) * prod(sqrt(b-a) ./ N);
%! psihat = sf2psihat(sf);
%! assert(psihatref(:),psihat(:))
