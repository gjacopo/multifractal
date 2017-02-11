%% ICIRCLE - Create circular structuring element
%
%% Syntax
%	S = ICIRCLE(R [, W])
%
%% Note
% Return a square matrix of zeros with a central circular region of radius R 
% of ones.  Matrix size is (2R+1) x (2R+1) or WxW.
% If R is a 2-element vector then it returns an annulus of ones, and the two 
% numbers are interpretted as inner and outer radii.

%% Function implementation
function s = icircle(r, w)

if ismatrix(r) 
  rmax = max(r(:));
  rmin = min(r(:));
else
  rmax = r;
end;

if nargin == 1,
  w = 2*rmax+1;
  c = rmax+1;
else
  c = ceil(w/2);
end;
s = zeros(w,w);

if ismatrix(r) 
  s = icircle(rmax,w) - icircle(rmin, w);
else
  [x,y] = find(s == 0);
  x = x - c;
  y = y - c;
  l = find(sqrt(x.^2+y.^2) <= r);
  s(l) = ones(length(l),1);
end;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function ret = ismatrix(M)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[sx sy] = size(M);
if sx>1 | sy>1
  ret=1;
else
  ret=0;
end;

