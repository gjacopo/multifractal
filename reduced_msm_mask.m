%% REDUCED_MSM_MASK - Compute the MSM and the chromatically reduced image.
%
%% Description
% Compute the MSM and the chromatically reduced image of a greyscale image 
% without taking into account pixels belonging to image mask.
%
%% Syntax
%    [MSM, dummy, Gx, Gy] = reduced_msm_mask(img, mask, ...
%                  upm_dens, upm_thres, [, flag]) 
%
%% Inputs
%    - img : original image,
%    - mask : image of mask (~=0 when belonging to the mask);
%      pixels of this mask are not taken into account in the
%      estimation of the MSM and in the reconstruction,
%    - upm_dens : density of MSM (0<upm_dens<1),
%    - upm_thres : threshold used to determine MSM,
%    - flag: if 1, display results,
%
%% Outputs
%    - dummy :chromatically reduced_from_msm image,
%    - MSM: most singular manifold, with gradient's orientation,
%    - Gx, Gy: derivatives of the reconstructed image.
%
%% Note
% The gradient that we retrieve here is the output of the following subsequent
% operations:
%    [gx,gy] = derive_spectral( img );
%    [gx,gy] = mask_gradient_msm( MSM, gx, gy,0 );
%    gx = propagation( gx, gy );
%    [Gx, Gy] = derive_spectral( gx );
%
%% See also
% Related:    
% reduced_msm 
% derive_spectral
% mask_gradient_msm
% propagation
% derive_spectral

%% Function implementation
function [MSM, dummy, Gx, Gy] = ...
    reduced_msm_mask(img, mask, upm_dens, upm_thres, flag)

if (exist('flag') ~= 1) flag =0; end;

[sx sy] = size(img);
[xeff yeff] = bits(sx,sy);

% 1) use of msm
% [MSM, gx, gy, dens] = msm( img,  gx, gy, ...
%                                upm_dens, upm_thres, flag );
% end of 1)

% replaced by: (avoid reccurent reprocessing)
% 2) use of derive_spectral and upm
 gx =zeros(xeff,yeff);
 gy =zeros(xeff,yeff);
 [ax, ay] = derive( img );
% le gradient des pixels du masque est mis a 0
  ax=ax.*(mask==0);
  ay=ay.*(mask==0); 
  gx(1:sx,1:sy)=ax;
  gy(1:sx,1:sy)=ay;
[gx, gy] = derive_spectral( img ); 

% de dimensions [xeff,yeff]=bits(sx,sy)
[MSM, dens] = upm_mask( img, mask, gx, gy, ...
				  upm_dens, upm_thres, flag );
modo=0;
% Computes unitary gradient: IFFT( FFT(\delta_{msm} * f#(-1)) ) 
[ax, ay] = derive_msm( MSM ); % de dimensions [xeff,yeff]
% Put the gradient to 0 for pixels outside the MSM and give a sign  
% to the MSM according to the orientation of the gradient
[gx, gy, MSM] = mask_gradient_msm( MSM, ax, ay, gx, gy, modo );
% same as mask_gradient for gx and gy
% end of 2)

% Operations pour le calcul de la PSNR 
% Reconstruction a partir du gradient sur la MSM: gx contient 
% l'image reconstruite 
gx = propagation( gx, gy );
PSNR = psnr( img, gx(1:sx,1:sy) );
fprintf('MSM with density %f at PSNR = %5f dB', dens, PSNR);

% Calcul du gradient de l'image reconstruite
[Gx, Gy] = derive_spectral( gx );
Gx = Gx(1:sx,1:sy); Gy = Gy(1:sx,1:sy);

% Calcul de l'image chromatiquement reduite
dummy =  dummy( MSM );
dummy = dummy(1:sx,1:sy);

if flag
  figure, imagesc(MSM), axis image, colormap gray;
  set(gca,'XTickLabel',{''}); set(gca,'YTickLabel',{''});
  title('Orientated Most Singular Manifold'), drawnow; 
  figure, imagesc(dummy), axis image, colormap gray;
  set(gca,'XTickLabel',{''}); set(gca,'YTickLabel',{''});
  title('Chromatically Reduced Image'), drawnow;
end;  
