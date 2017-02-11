%% MASK_GRADIENT - Apply the mask of MSM on the gradient images. 
%
%% Syntax
%       [Gx,Gy] = mask_gradient( msm, gx, gy )
%
%% Note
% gx=gy=0 on pixels outside the MSM (ie. when, MSM=C0) and outside the dimension 
% of the MSM.
%
%% See also
% Related:    
% derive_msm_unitary
% upm

%% Function implementation
function [Gx,Gy] = mask_gradient( msm, gx, gy )

C0=255; CP=0; CM=127;
% C0 = grey level associated to 0
% CM = grey level associated to -1  
% CP = grey level associated to 1

[sx sy] = size(msm);
[xeff yeff] = size(gx);
% we must have: [xeff, yeff] = bits(sx,sy);

Gx = zeros(xeff, yeff);
Gy = zeros(xeff, yeff);
Gx(1:sx,1:sy) = (msm ~= C0) .* gx(1:sx,1:sy); % +(msm==C0)*0.
Gy(1:sx,1:sy) = (msm ~= C0) .* gy(1:sx,1:sy); % +(msm==C0)*0.
