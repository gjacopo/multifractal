%% VEC_NORMALISE - Pixelwise normalisation of an image.
%
%% Description
% Normalise a vectorial field in the complex plane.
%
%% Syntax
%        [ax, ay]= vec_normalise(ax, ay)
% 
%% See also
% Related:
% <|derive_msm|>

%% Function implementation
function [ax, ay]= vec_normalise(ax, ay)

Mod = sqrt( ax.*ax + ay.*ay);
Mod = (Mod>1e-30) .* Mod + (Mod<=1e-30);
ax = ax ./ Mod;
ay = ay ./ Mod;
  
