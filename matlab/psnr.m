%% PSNR - Compute the PSNR ratio between an image and its approximation.
%
%% Syntax
%        PSNR = psnr(img, err)
%
%% Note
% Signal-to-noise (SNR) measures estimates the quality of the reconstructed image 
% compared with the original image.
% Computation of the mean squared error MSE of the reconstructed image IR as 
% follows:
%          MSE = \frac{\sum_x (IR(x) - I(x))^2} {NM}
% where I is the original image and (N,M) is the size of the image
% 1) Compute the difference between original and reconstructed images
% 2) Compute the root mean squared error RMSE which is the square root of MSE:
%          RMSE = \sqrt{MSE} 
% 3) Compute the PSNR as follows:
%          PSNR = 20 \log_10 (\frac{255.}{RMSE})
% (or similarly: PSNR = - 10 \log_10 (\frac{MSE}{255.^2}), 
% depending on definitions).
% Note that error metrics are computed on the luminance signal only so the pixel 
% values I and IR range between black (0) and white (255).
%
%% See also
% Related:    
% msm 
% reduced_msm

%% Function implementation
function PSNR = psnr(img, err)

errmin=min(min(err));
errmax=max(max(err));
wds=errmax-errmin;
if (wds>1e-30)
  eerr = 255.*(err-errmin) / wds;   %  eerr = err;
end;
immin=min(min(img));
immax=max(max(img));
wds=immax-immin;
if (wds>1e-30)
  iim = 255.*(img-immin) / wds;  % iim =img;
end;

eerr = iim - eerr;

% Normalisation du signal dans [0,255]
s =length(eerr(:));
% mse = sum(eerr(:).^2)
rmse = sqrt(sum(eerr(:).^2) / s);
PSNR = 20. * log(255./rmse) / log(10.);

% Normalisation du signal dans [0,255]
%  img=normalizeMat(img,0,255);
%   err=normalizeMat(err,0,255);
%   s =length(img(:));
%   mse = sum(sum((err-img).^2))
%   rmse = sqrt(sum(sum((err-img).^2)) / s)
%   PSNR = 20. * log(255./rmse) / log(10.)


