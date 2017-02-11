%% DERIVE_NORM - Compute the norm of a naive derivative.
%
%% Syntax
%        der = derive_norm( img )

%% Function implementation
function der = derive_norm( img )

[sy sx]=size(img);
img_trans = [ img, img(:,1); img(1,:), 0.];

dX = img_trans(1:sy,2:sx+1) - img;
dY = img_trans(2:sy+1,1:sx) - img;

der = sqrt(dX.*dX + dY.*dY);
