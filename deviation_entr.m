function dev = deviation_entr(im,entr)

%%% 
[gimy, gimx] = derive( im );
[gentry, gentrx] = derive( entr );
% [gimy, gimx] = derive_spectral( im );
% [gentry, gentrx] = derive_spectral( entr );
dev=asigna_orientation(gimy, gimx, gentry, gentrx);
figure,imagesc(dev),colormap jet, axis image

return;
[sx,sy]=size(im);
[gx,gy] = derive_spectral(entr);

%%% 
% UNITARY - Computes the MSM and the chromatically 
% reduced_from_msm image of an image with a naive unitary gradient instead 
% of the true one (but with the same direction and orientation).
upm_dens=0; upm_thres=1.; 
flag=0;
[MSM,dummy] = unitary(im, upm_dens, upm_thres, flag);

d = calcula_deviation (dummy,gx,gy);
d = d(1:sx,1:sy);
figure,imagesc(d),colormap jet, colorbar,axis image

segdev= 0*(d<pi/8) + (pi/8)* ((d>=pi/8) & (dev<pi/4)) ...
	+ (pi/4)* ((d>=pi/4) & (d<3*pi/8)) ...
	+ (3*pi/8)* ((d>=3*pi/8) & (d<pi/2)) ...
	+ (pi/2)* ((d>=pi/2) & (d<5*pi/8)) ...
	+ (5*pi/8)* ((d>=5*pi/8) & (d<3*pi/4)) ...
	+ (3*pi/4)* ((d>=3*pi/4) & (d<7*pi/8)) ...
	+ (7*pi/4)* ((d>=7*pi/4) & (d<pi));
	
% segdev= (pi/8)* ((d<pi/8) | (d<7*pi/8)) ...
% 	+ (pi/4)* ( ((d>=pi/8)&(d<pi/4)) | ((d>=3*pi/4)&(d<7*pi/8)) ) ...
% 	+ (3*pi/4)* ( ((d>=pi/4)&(d<3*pi/8)) | ((d>=5*pi/8)&(d<3*pi/4)) ) ...
% 	+ (pi/2)* ((d>=3*pi/8) & (d<5*pi/8)) ;
segdev= 0.8* ((d<pi/8) | (d<7*pi/8)) ...
	+ 1.4* ( ((d>=pi/8)&(d<pi/4)) | ((d>=3*pi/4)&(d<7*pi/8)) ) ...
	+ 1.7 * ( ((d>=pi/4)&(d<3*pi/8)) | ((d>=5*pi/8)&(d<3*pi/4)) ) ...
	+ 2.5 * ((d>=3*pi/8) & (d<5*pi/8)) ;
	
segdev=  ((d>=3*pi/8) & (d<5*pi/8)) & (im<220);


function dev = deviation (dummy,gx,gy)

EXP_MU=1.;
gx = filter_spectral(gx,0.,-EXP_MU);
gy = filter_spectral(gy,0.,-EXP_MU);

ax = filter_spectral(dummy,0.,-EXP_MU);
[ax,ay] = derive_spectral(ax);
dev = orientation(ax,ay,gx,gy);
  
function orient=orientation(ax,ay,gx,gy)

moda=sqrt(ax.*ax+ay.*ay);
modg=sqrt(gx.*gx+gy.*gy);
mod=modg.* moda;
% mod = mod .* (mod~=0) + (mod==0);

prod=gx.*ax + gy.*ay;
prod = prod ./ mod;

orient=acos(prod);
