%% Copyright (C) 2016 Marco Caliari
%% Copyright (C) 2010 VZLU Prague, a.s., Czech Republic
%%
%% Author: Jaroslav Hajek
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; see the file COPYING.  If not, see
%% <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn{Function File} {@var{y} =} ndcovlt (@var{x}, @var{t1}, @var{t2}, @dots{})
%% @deftypefnx{Function File} {@var{y} =} ndcovlt (@var{x}, @var{t})
%% Computes an n-dimensional covariant linear transform of an n-d tensor, given a
%% transformation matrix for each dimension. The number of columns of each transformation
%% matrix must match the corresponding extent of @var{x}, and the number of rows determines
%% the corresponding extent of @var{y}. For example:
%%
%% @example
%%   size (@var{x}, 2) == columns (@var{t2})
%%   size (@var{y}, 2) == rows (@var{t2})
%% @end example
%%
%% The element @code{@var{y}(i1, i2, @dots{})} is defined as a sum of
%%
%% @example
%%   @var{x}(j1, j2, @dots{}) * @var{t1}(i1, j1) * @var{t2}(i2, j2) * @dots{}
%% @end example
%%
%% over all j1, j2, @dots{}. For two dimensions, this reduces to
%% @example
%%   @var{y} = @var{t1} * @var{x} * @var{t2}.'
%% @end example
%%
%% [] passed as a transformation matrix is converted to identity matrix for
%% the corresponding dimension. A single cell @var{t} containing the
%% transformation matrices can be used.
%%
%% @end deftypefn

function y = ndcovlt (x, varargin)
% Function file: y = ndcovlt(x,varargin)
%
% Description: Compute an n-dimensional linear transform of an nd tensor.
  nd = max (ndims (x), nargin - 1);

  if (iscell (varargin{1}))
     varargin = varargin{1};
  end

  % check dimensions
  for i = 1:nd
    ti = varargin{i};
    if (isnumeric (ti) && ndims (ti) == 2)
      [r, c] = size (ti);
      if (r + c == 0)
        varargin{i} = eye (size (x, i));
      elseif (c ~= size (x, i))
        error ('ndcovlt: dimension mismatch for x-th transformation matrix');
      end
    else
      error ('ndcovlt: transformation matrices must be numeric 2d matrices');
    end
  end

  ldp = [2:nd, 1];
  %% First transformation.
  y = ldtrans (x, varargin{1});

  %% Always shift one dimension.
  for i = 2:nd-1
    y = ldtrans (permute (y, ldp), varargin{i});
  end

  %% Permute to normal order now to save one permutation.
  if (nd > 2)
    y = ipermute (y, [nd-1:nd, 1:nd-2]);
  end

  %% Now multiply from the right.
  szy = size (y);
  szy(end+1:nd-1) = 1;
  m = varargin{nd};
  szy(nd) = size (m, 1);
  y = reshape (y, [], size (y, nd));
  y = reshape (y * m.', szy);

end

function y = ldtrans (x, m)
  sz = size (x);
  sz(1) = size (m, 1);
  y = reshape (m * x(:,:), sz);
end
