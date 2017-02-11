%% UPM - Compute the MSM and the density of the MSM.
%
%% Description
% Compute the MSM and the density of the MSM in the image according to 
% some adjusting parameters.
%
%% Syntax
%  [msm, dens] = upm( img, gx, gy, upm_dens, upm_thres [, flag] )
% 
%% Inputs
%      - img: original image,
%      - gx, gy: derivatives of img in x- and y-directions
%        (calculated with derive_spectral: 
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
    upm( img, gx, gy, upm_dens, upm_thres, flag  )

if (exist('flag') ~= 1) flag =0; end;
 
[sx sy] = size(img); % sx: nb lignes, sy: nb colonnes

disp = distribution_upm( img, gx, gy );
% DEBUG: disp=img;

media = sum(sum(disp)) / (sx*sy);
% Calcul du seuil de determination des valeurs des exposants de la MSM 
% a partir de leur distribution 
if upm_dens>0.
  [umbral, n] = quantile_threshold( disp, upm_dens, flag );
  umbral = umbral / media;
else
  umbral=upm_thres;
  Nhisto=10000;
  n = hist(disp(:),Nhisto) / (sx*sy);
end;

[msm, dens] = msm_from_sing( disp, media, umbral );

% done upper in the calling functions:
% [Gx, Gy] = mask_gradient( msm, gx, gy );
% err = propagation( Gx, Gy );
% PSNR = psnr( img, err(1:sx,1:sy) );
% fprintf('MSM with density %f at PSNR = %5f dB', dens, PSNR);

