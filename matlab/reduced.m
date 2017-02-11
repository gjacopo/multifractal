%% REDUCED_FROM_MSM - Perform the reconstruction of the image. 
%
%% Description
% Perform the reconstruction of the image using the reduced unitary essential 
% vectorial field, which :
%   - is unitary on the MSM, null otherwise
%   - is perpendicular to the MSM on each pixel of the MSM,
%   - has identical orientation as the original gradient.
%
%% Syntax
%       dummy =  dummy( MSM )
% 
% Note 
% The orientation of the gradient is given in the msm.
%
%% See also
% Related:    
% reduced_msm 

%% Function implementation
function dummy =  dummy( MSM )

% Computes the reduced_from_msm unitary essential vectorial field
[Gx, Gy] = derive_msm_unitary( MSM );

% DEBUG
% [Gx1, Gy] = derive_msm_unitary( MSM );
% [Gx,Gy] = derive_spectral( Gx1 );
% [Gx,Gy] = mask_gradient( MSM, Gx, Gy );
% ENDDEBUG

% Reconstructs image
dummy = propagation(Gx,Gy);

% DEBUG
% subplot(1,3,2), imshow(dummy,[]);
% subplot(1,3,3), imshow(Gx1-dummy,[]);
% ENDDEBUG
