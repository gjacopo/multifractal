% SOURCE_FROM_MSM - Compute the chromatically reduced images and the source. 
%
%% Description
% Compute the chromatically reduced images and the source from an image and 
% its MSM. Use parts of code lines include in file reduced_MSM. Avoid 
% re-processing for computing the MSM.
%
%% Syntax
%    [Mx,My,MM,orMM,dummy] = source_from_msm( img, MSM, flag )

%% Function implementation
function [MM, orMM, Mx, My, dummy] = source_from_msm( img, MSM, flag )

C0=255; CP=0; CM=127;
% C0 = grey level associated to 0
% CM = grey level associated to -1
% CP = grey level associated to 1

[sx sy] = size( MSM );
[xeff, yeff] = bits(sx,sy);

% Computes the reduced_from_MSM unitary essential vectorial field
% and then computes the chromatically reduced_from_MSM image.
dummy =  dummy( MSM );    
dummy = dummy(1:sx,1:sy);

% if flag
%   figure, dispImg(MSM), colormap gray, 
%   title('Orientated Most Singular Manifold');
%   figure, imagesc(dummy), colormap gray, 
%   title('Chromatic image');
% end;

% On recupere le gradient comme en sortie de reduced_MSM
[gx, gy] = derive_spectral( img );
[gx,gy] = mask_gradient( MSM, gx, gy );
gx = propagation( gx, gy );
[Gx, Gy] = derive_spectral( gx );
gx = Gx(1:sx,1:sy); 
gy = Gy(1:sx,1:sy);

silog=1;
[MM, orMM, Mx, My] = ...
    source( dummy, gx, gy, silog, flag );

% [ax, ay] = derive_spectral( dummy );
% ax = ax(1:sx,1:sy); 
% ay = ay(1:sx,1:sy);
% aMod = sqrt(ax.^2 + ay.^2);
% gMod = sqrt(gx.^2 + gy.^2);
% angl = (ax.*gx + ay.*gy)  ./ (aMod .* gMod);
% figure, imagesc(angl), colormap jet;

% Mx = 0;
% My = 0;
% orMM = angl;
