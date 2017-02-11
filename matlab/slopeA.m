%% SLOPE_SPECTRAL - Compute the \eta coefficient so that A(f) varies
% like f^{\eta}.
%
%% Syntax
%      slope = slope_spectral(A, nbins) 
%
%% Description
% Check that A^2(f) varies in f^\eta since it is defined in frequency domain:
%          A(f) = |\vec{v}_\infty \cdot \vec{f}| / f
% where \vec{f} is the frequency vector, f its module, and \vec{v}_\infty is 
% the essential gradient and \cdot is the scalar product.
% As expressed in equation (21) of:
%   "Reconstructing images from their most singular fractal manifold"
%   Turiel et del Pozo
% the following holds:
%          S(\vec{f}) = g^2(f) . A^2(vec{f})

%% Function implementation
function slope_spectral(A, nbins) 

% if (exist('freq')~=1)  freq=; end;

[xeff yeff] = size(A);
if (exist('nbins') ~= 1)  nbins = 400;  end;

[X Y] = meshgrid( [0:yeff/2,0,-yeff/2+1:-1]/yeff, ...
		  [0:xeff/2,0,-xeff/2+1:-1]/xeff);
% Vecteur frequence
dX = sin(pi.*X);
dY = sin(pi.*Y);
% Vecteur frequence
freq = sqrt(dX.^2 + dY.^2); 
size(freq)
% or:  dX =X; dY =Y;
% freq = sqrt(X.^2 + Y.^2);
freq = round(freq/max(max(freq))*(nbins-1)) + 1;


amp = zeros(1, nbins);
fcount = ones(1, nbins);
for r= 1:xeff
  for c = 1:yeff
    ind = freq(r,c);
    amp(ind) = amp(ind) + A(r,c);
    fcount(ind) = fcount(ind)+1;
  end;
end;

amp = amp ./ fcount;             
f = [1:nbins] / nbins * sqrt(.5); 

cut=2;
verb=1;
slope = specSlope( amp, f, cut, verb);
