%% REDUCED_UNITARY - Compute a chromatically reduced image.
%
%% Description
% Compute the chromatically reduced image reconstructed using an unitary
% gradient instead of the original one (but with the same direction and 
% orientation).
%
%% Syntax
%      dummy = reduced_unitary(gx,gy)
% 
%% See also
% Related:     
% unitary

%% Function implementation
function dummy = reduced_unitary(gx,gy)

% Calcul du  module du gradient
mod = sqrt(gx.^2 + gy.^2);
mod = (mod>1.e-30) .* mod + (mod<=1.e-30);

ax = gx ./ mod;
ay = gy ./ mod;

dummy = propagation(ax,ay);
