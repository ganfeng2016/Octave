## Copyright (C) 2015, Feng Gan <cesgf@mail.sysu.edu.cn;sysucesgf@163.com>
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
## 
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.
##

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{C}] =} helps (@var{X}, pcn, selInfo)
## Window factor analysis with the model @code{X = CS}.
##
## @itemize
## @item
## @code{X} is a spectral matrix whose each column is a spectrum.
## @item
## @code{pcn} is the number of components.
## @item
## @code{selInfo} is a matrix containing selective regions and zero concentration regions.
## @end itemize
##
## Return values
##
## @itemize
## @item
## @code{C} is a concentration matrix
## @end itemize
## @end deftypefn
##
## References:
##
## - Malinowski, E. R. Window factor analysis: Theoretical derivation and application to flow injection analysis data. J. Chemom. 6, 29–40 (1992).
## - Den, W. E. & Malinowski, E. R. INVESTIGATION OF COPPER ( I1 ) -ETHYLENEDIAMINETETRAACETATE COMPLEXATION BY WINDOW FACTOR ANALYSIS OF ULTRAVIOLET SPECTRA *. 7, 89–98 (1993).
##
## Author:  Feng Gan
## Latest revision: 2015-10-17
## Create: 1999-01-20

function [C] = wfa(X, pcn, selInfo)

if nargin < 3
	error('You need three input parameters!');
end

C = [];

[m,n] = size(selInfo);
for i = 1:m

	X1 = X(:, selInfo(i,1):selInfo(i,2));
	X2 = X(:, selInfo(i,3):selInfo(i,4));
	D = [X1 X2];

	[S0, L0, C0] = svd(D);
	S0 = S0(:,1:pcn-1);
	V = S0 * S0';
	I = eye(size(V));
	M = (I - V) * X;

	[u,s,v] = svd(M');
	ci = u(:,1);
	[maxv,locat] = max(abs(ci));
	if ci(locat) < 0
		ci = -1.0 * ci;
	end
    C = [C ci];

endfor

endfunction

%!demo
%! load ./Data/masartdata.mat
%! X = X';
%! selInfo = [1 8 23 50];
%! pcn = 5;
%! [C] = wfa(X, pcn, selInfo);
%! plot(C)

