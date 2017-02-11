%% VEC_DIVIDE - Pixelwise division of an image.
%
%% Description
% Compute the complex division of two vectorial fields seen as components 
% in the complex plane.
%
%% Syntax 
%     [Gx, Gy] = vec_divide( ax, ay, gx, gy)
%
%% Outputs
%     - Gx = (ax.*gx + ay.*gy) / mod, 
%     - Gy = (ax.*gy - ay.*gx) / mod,
% where mod stands for the norm of (ax,ay).
%
%% See also  
% Related:
% <|source|>

%% Function implementation
function [Gx, Gy] = vec_divide( ax, ay, gx, gy) 

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Remarque:
% Si f=(Gx,Gy) est vu comme une fonction de valeurs complexes,
% alors, en supposant f R-derivable, f est holomorphe, puisque
% ...
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Mod = ax.^2 + ay.^2;
Mod = (Mod>1e-30) .* Mod  + (Mod<=1e-30).*1;

% Division complexe: z1 / z2 = (z1*\bar{z2}) / (z2*\bar{z2})
%                            = (z1*\bar{z2}) / |z2|^2
Gx = (ax.*gx + ay.*gy) ./ Mod;
Gy = (ax.*gy - ay.*gx) ./ Mod;


