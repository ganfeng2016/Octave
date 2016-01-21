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
## @deftypefn {Function File} {[@var{s1}, @var{s2}, @var{prl}, @var{isp}, @var{p}] =} smcr (@var{X})
## Self-modeling curve resolution for only two component system.
## with the model @code{X = CS}.
##
## @itemize
## @item
## @code{X} is a spectral matrix whose each row is a spectrum.
## @end itemize
##
## Return values
##
## @itemize
## @item
## @code{s1} is a spectral matrix of the first component.
## @item
## @code{s2} is a spectral matrix of the first component.
## @item
## @code{prl} is a matrix that contains the potential regions.
## @item
## @code{isp} is a matrix that contains the intersection points.
## @item
## @code{p} is a vector that contains coefficients.
## @end itemize
## @end deftypefn
##
## References:
## - W.E. Lawton, E.A. Sylvester, Self-modeling curve resolution, Technometrics. 13 (1971) 617 - 633.
## - 梁逸曾，俞汝勤，化学计量学，高等教育出版社：北京，2003
##
## Lastest Revision: 2015-07-21
## Create date: 1999-01-16

function [s1,s2,prl,isp,p] = smcr(X)

s1 = [];
s2 = [];
[m,n] = size(X);
for i = 1:m
   X(i,:) = X(i,:)/sum(X(i,:));
endfor
[u,s,v] = svd(X);
if max(v(:,1)) < 0
    u = -1.0.*u;
    v = -1.0.*v;
endif
u = u * s;
max_u1 = max(u(:,1));
vr = v(:,1) ./ v(:,2);
ur = u(:,2) ./ u(:,1);
[mvr,nvr] = size(vr);
vrz = [];
vrf = [];
for i = 1:mvr
   if v(i,2) > 0
      vrz = [vrz;vr(i)];
   else
      vrf = [vrf;vr(i)];
   endif
endfor

## calculate slopes
slop_L1 = - min(vrz);
slop_H1 = min(abs(vrf));
slop_L2 = min(ur);
slop_H2 = max(ur);

## calculate potential regions
t1 = 0.0:0.01:(max_u1 + 0.2);
L1 = slop_L1 * t1;
H1 = slop_H1 * t1;
L2 = slop_L2 * t1;
H2 = slop_H2 * t1;
prl = [t1; L1; H1; L2; H2];

## caculate intersection points
b = sum(v(:,1));
a = sum(v(:,2));
p = [a b];
[L1zb] = fc(slop_L1,a,b);
[H1zb] = fc(slop_H1, a, b);
[L2zb] = fc(slop_L2,a,b);
[H2zb] = fc(slop_H2, a, b);
isp = [L1zb;H1zb;L2zb;H2zb];

## calculate spectra
bp = 0.0;
ep = 0.0;
if H1zb(1) < H2zb(1)
   bp = H1zb(1);
   ep = H2zb(1);
else
   bp = H2zb(1);
   ep = H1zb(1);
endif
for i = bp:0.001:ep
   t2 = (1 - b * i) / a;
   s1tmp = i * v(:,1) + t2 * v(:,2);
   s1 = [s1 s1tmp];
endfor
if L1zb(1) < L2zb(1)
   bp = L1zb(1);
   ep = L2zb(1);
else
   bp = L2zb(1);
   ep = L1zb(1);
endif
bc = 0.0;
if (a/b) < 0
   p1 = ep;
   p2 = bp;
   bc = -0.0001;
else
   p1 = bp;
   p2 = ep;
   bc = 0.0001;
endif
for i = p1:bc:p2
   t2 = (1 - b * i) / a;
   s2tmp = i*v(:,1) + t2*v(:,2);
   s2 = [s2 s2tmp];
endfor

endfunction

##
function [zb] = fc(k, a, b)
x = 1/(a*k + b);
y = k/(a*k + b);
zb=[x y];
endfunction

%!demo
%! load ./Data/smcrm.mat
%! [s1,s2,prl,isp,p] = smcr(X);
%! plot(s1)
%! hold on
%! plot(s2)


