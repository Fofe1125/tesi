function sf = sfrun(pdb,eta,ic,T,nsteps,nfiles,rfn,countr)
% Function file: sf = sfrun(pdb,eta,ic,mirror,T,nsteps,nfiles,rfn)
% Function file: sf = sfrun(pdb,eta,ic,mirror,T,nsteps,nfiles,rfn,countr)
%
% Description: Run a complete superfluid simulation.
%
% Input variables:
%   pdb:    physical domain boundaries in the axis form
%   eta:    row vector of numbers of grid points in the physical domain,
%           eta(d) even
%   ic:     function handle defining the structure (see below) of the initial
%           condition, in the form sf = ic(x), where x is the cell
%           array containing the vectors x{d}
%   T:      final integration time
%   nsteps: number of time steps
%   nfiles: number of output mat-files, in addition to the first file,
%           containing the structure (see below) of the solutions at
%           initial time t0, t0+(T/nsteps)*(nsteps/nfiles),
%           t0+2*(T/nsteps)*(nsteps/nfiles), ...,
%           t0+nfiles*(T/nsteps)*(nsteps/nfiles)
%   rfn:    root of output file name (rfn000.mat, rfn001.mat, ...),
%           where the number of leading zeros is automatically
%           computed by nfiles and the optional input countr
%   countr: optional starting value for file names (default is 0)
%
% Output variables:
%   sf:     structure defining the function psi at time T, where
%           sf.pdb:    are the physical domain boundaries
%           sf.psipdb: is the complex nd-array of psi in the physical domain
%           sf.mirror: is the row vector of mirroring flags true/false
%           sf.t:      is the simulation time
dd = length(eta);
x = mkgrid(pdb,eta); % make grid in the physical domain
sf = ic(x); % initial condition
t = sf.t;
if (nargin == 7)
  countr = 0;
end
[psi,a,b,Lambda2] = sfpregpe(sf); % preproc.
[mass,energy] = sfEm(sf,Lambda2); % mass and energy
fprintf(1,'Initial mass:   %.6e\n',mass);
fprintf(1,'Initial energy: %.6e\n',energy);
basefilename = sprintf('%s%%0.%dd.mat',rfn,floor(log10(nfiles+countr)+1));
filename = sprintf(basefilename,countr);
save('-v6',filename,'sf')
fprintf(1,'sf at time %g written on %s\n',t,filename);
k = T/nsteps; % time step
jout = nsteps/nfiles; % number of time steps between output files
dT = k*jout; % integration time between output files
idx = cell(1,dd);
for d = 1:dd
  idx{d} = 1:eta(d);
end
for i = 1:nsteps/jout
  psi =  sfgpe(psi,dT,jout,Lambda2); % dT integration time
  t = t+dT;
  sf.psipdb = psi(idx{:}); % extract solution on pdb
  sf.t = t;
  countr = countr+1;
  filename = sprintf(basefilename,countr);
  save('-mat',filename,'sf')
  fprintf(1,'sf at time %g written on %s\n',t,filename);
end
[mass,energy] = sfEm(sf,Lambda2);
fprintf(1,'Final mass:     %.6e\n',mass);
fprintf(1,'Final energy:   %.6e\n',energy);
