% Description: Script for selecting NFFT usage.

% This script automatically recognizes whether Matlab or GNU Octave
% is running, however the user has to set the following two variables.

% Set the path of the NFFT Matlab interface (if available),
% otherwise set matlabpath = '';

matlabpath = '/usr/local/src/nfft-3.3.2-matlab/matlab/nfft';

% Set the path of the NFFT GNU Octave interface (if available),
% otherwise set octavepath = '';

octavepath = '/usr/local/src/nfft-3.3.2-octave/matlab/nfft';

% If NFFT is not available (deprecated) or the user does not want to
% use it (e.g., testing purposes), all paths must be set to ''.

if (exist('OCTAVE_VERSION'))
  if (isempty(octavepath))
    havenfft = false;
  else
    havenfft = true;
    addpath(octavepath)
  end
else
  if (isempty(matlabpath))
    havenfft = false;
  else
    havenfft = true;
    addpath(matlabpath)
  end
end
%!demo
%! nfftpath
%! simple_test
