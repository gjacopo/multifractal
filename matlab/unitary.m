%% UNITARY - Compute the MSM and the chromatically reduced image of an image.
%
%% Description
% Compute the MSM and the chromatically reduced image of an image with 
% a naive unitary gradient instead of the true one (but with the same 
% direction and orientation).
%
%% Syntax 
%   [MSM,dummy] = unitary(img, upm_dens, upm_thres, flag)
% 
%% Inputs
%    - img : original image,
%    - upm_dens : density of MSM (0<upm_dens<1),
%    - upm_thres : threshold used to determine MSM,
%    - flag: if 0, display results,
%    - MSM: most singular manifold,
%    - dummy :chromatically reduced image with naive gradient,
%      ie. unitary and with the same direction as the gradient
%      of the image.
%
%% Note:
% Dimension of dummy are powers of 2.)

%% Function implementation
function [MSM,dummy] = unitary(img, upm_dens, upm_thres, flag)

[sx sy] = size(img);
[gx,gy] = derive_spectral( img );

dummy = reduced_unitary(gx,gy);
d=dummy(1:sx,1:sy);

[MSM, Gx, Gy, dens] = msm( d, upm_dens, upm_thres );

err = propagation( Gx, Gy );
PSNR = psnr( img, err(1:sx,1:sy) );
fprintf('MSM with density %f at PSNR = %5f dB', dens, PSNR);

if flag
  figure, colormap gray, imagesc(MSM), axis image, 
  title('Orientated MSM,  unitary'), drawnow;
  figure, colormap gray, imagesc(d), axis image, 
  title('Chromatically reduced Image, unitary'), drawnow;
end;
