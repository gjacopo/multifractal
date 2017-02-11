%% UPM_MASK - Compute the MSM and the density of the MSM.
%
%% Description
% Compute the MSM and the density of the MSM in the image taking into account a 
% mask for which pixels have not to be considered and according to some adjusting 
% parameters.
%
%% Syntax
%  [msm, dens] = upm_mask( img, mask, gx, gy, upm_dens, upm_thres [, flag] )
%
%% Inputs
%      - img: original image,
%      - mask: mask of pixels which should not be taken into
%        account in the  computation of the MSM; in any case,
%        such pixels cannot belong to the MSM (see impainting)
%        (~=0 when belonging to the mask),
%      - gx, gy: derivatives of img in x- and y-directions
%        (calculated with derive_spectral.
%
%% Outputs
%      - msm: mask of the MSM,
%      - dens: density of the MSM in the image.
% 
%% See also
% Related:     
% reduced_msm
% msm

%% Function implementation
function [msm, dens] = ...
    upm_mask( img, mask, gx, gy, ...
			upm_dens, upm_thres, flag  )




if (exist('flag') ~= 1) flag =0; end;
 
[sx sy] = size(img); % sx: nb lignes, sy: nb colonnes

disp = distribution_upm( img, gx, gy );

% the singularity of pixels belonging to the mask (whom value is ~=0
% in the mask) are set to 0.
disp = disp .* (mask==0);

% number of pixels in the mask
nbm = size(find(mask~=0),1);

media = sum(sum(disp)) / (sx*sy-nbm);
% Calcul du seuil de determination des valeurs des exposants de la MSM 
% a partir de leur distribution 
if upm_dens>0.
  [umbral, n] = quantile_threshold( disp, upm_dens, flag );
  umbral = umbral / media;
else
  umbral=upm_thres;
  Nhisto=10000;
  T = disp(:);
  n = hist(T,Nhisto) / (sx*sy-nbm);
end;

% Pixels of the mask are then set to the min of disp so that they
% cannot be detected as belonging to the MSM
%
disp = disp .* (mask==0) + min(min(disp)) .* (mask~=0);
[msm, dens] = msm_from_sing( disp, media, umbral );

% done upper in the calling functions:
% [Gx, Gy] = mask_gradient( msm, gx, gy );
% err = propagation( Gx, Gy );
% PSNR = psnr( img, err(1:sx,1:sy) );
% fprintf('MSM with density %f at PSNR = %5f dB', dens, PSNR);

