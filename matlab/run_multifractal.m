% RUN_MULTIFRACTAL - Launch multifractal decomposition, source extraction 
% and reconstruction.
%
%% Syntax
%   [MSM, dummy, source, Gx, Gy, signal] = ...
%               run_multifractal(flag [, img] )

%% Function implementation
function [MSM, dummy, source, Gx, Gy, signal] = ...
    run_multifractal(flag, img)

if exist('flag') ~= 1 flag=0; end;

upm_dens=0.45; upm_thres=1.;
[img,med] = anorma(img);
[sx,sy]=size(img);

% 1) First computes the MSM and the chromatically reduced_from_msm image of the
% image with a naive unitary gradient insted of the true one (but with 
% the same direction and orientation).
if flag 
  fprintf('\n1) Computes the unitary MSM and the unitary chromatically reduced_from_msm image');
end;
[unitaryMSM,unitaryMSMdummy] = ...
    unitary(img, upm_dens, upm_thres, flag);
unitaryMSMdummy=unitaryMSMdummy(1:sx,1:sy);

% 2) Then computes the MSM and the chromatically reduced_from_msm image of the
% image (associated with the true gradient of image).
if flag 
  fprintf('\n2) Computes the orientated MSM and the chromatically reduced_from_msm image');
end;
[MSM, dummy, Gx, Gy] = ...
    reduced_msm(img, upm_dens, upm_thres, flag);

% 3) Computes the source of an image, starting from the estimation
% of the MSM and of the derivatives computed by reduced_msm.
silog=1;
if flag 
  fprintf('\n3) Computes the source associated to the MSM');
end;
[source, source_vec, Gx, Gy] = ...
    source( dummy, Gx, Gy, silog, flag );

% 4) Computes the reconstruction of the image.
if flag 
  fprintf('\n4) Computes the reconstruction of the image');
end;
signal = reconstruction(dummy, Gx, Gy, flag );
signal = shift(signal,med);
