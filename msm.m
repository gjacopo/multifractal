%% MSM - Compute the MSM.
%
%% Syntax
%   [MSM, Gx, Gy, dens] = msm( img, upm_dens, upm_thres[, flag] )
%
%% Inputs
%    - img : original image,
%    - upm_dens : density of MSM (0<upm_dens<1),
%    - upm_thres : threshold used to determine MSM,
%    - flag: if 0, display results.
%
%% Outputs
%    - MSM: most singular manifold, (without information about 
%      the gradient's orientation),
%    - Gx, Gy: derivatives of the original image on which the
%      mask of MSM is applied.
%
%% See also
% Related:    
% reduced_msm 
% unitary

%% Function implementation
function [MSM, Gx, Gy, dens] =  msm( img, upm_dens, upm_thres, flag )

if (exist('flag') ~= 1) flag =0; end;
 
[gx, gy] = derive_spectral( img ); % de dimensions [xeff,yeff]=bits(sx,sy)

[MSM, dens] = upm( img, gx, gy, ...
			     upm_dens, upm_thres, flag );

% done elsewhere
% [sx sy] = size(img);
% [Gx, Gy] = mask_gradient( MSM, gx, gy );
% err = propagation( Gx, Gy );
% PSNR = psnr( img, err(1:sx,1:sy) );
% fprintf('MSM with density %f at PSNR = %5f dB', dens, PSNR);

modo=0;
% Computes unitary gradient: IFFT( FFT(\delta_{msm} * f#(-1)) ) 
[ax, ay] = derive_msm( MSM ); % de dimensions [xeff,yeff]
% Put the gradient to 0 for pixels outside the MSM and give a sign  
% to the MSM according to the orientation of the gradient
[Gx, Gy, MSM] = mask_gradient_msm( MSM, ax, ay, gx, gy, modo );
