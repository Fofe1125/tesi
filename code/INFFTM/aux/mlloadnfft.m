% Description: Script to fix bug. no.961694 in certain Matlab versions.
% In order to prevent Matlab bug no.961694, nfftmex has to be loaded
% as soon as possible. It can be done by running the following script
% as soon as Matlab starts.
nfftpath
plan = nfft_init_3d(10,10,10,10);
nfft_finalize(plan);
