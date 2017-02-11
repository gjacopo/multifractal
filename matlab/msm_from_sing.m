%% MSM_FROM_SING - Compute the MSM and its density according to the
% distribution of the singularity exponents.
%
%% Syntax:
%      [msm, dens] = msm_from_sing( disp, media_tot, umbral )
%
%% See also
% Related:    
% upm

%% Function implementation
function [msm, dens] = msm_from_sing( disp, media_tot, umbral )

C0=255.; CP=0.; 
% C0 = grey level associated to 0
% CP = grey level associated to 1  

[sx sy] = size(disp);

% msm = ones(sx,sy) * C0;
msm = (disp >=umbral*media_tot) .* CP ...
      + (disp <umbral*media_tot) .* C0;

dens = size(find(msm~=C0),1) / (sx*sy);
% ou dens=1-sum(sum(msm==C0))/(sx*sy*C0), mais plus dangereux...
