%% RECONSTRUCTION - Compute the reconstruction of an image.
%
% Syntax
%        signal = reconstruction(dummy, gx, gy [, flag] )
%
%% Inputs
%   - dummy :chromatically reduced_from_msm image (computed by 
%     reduced_msm or unitary),
%   - gx, gy: images of the component of the source computed 
%     by source.

%% Function implementation
function signal = reconstruction(dummy, gx, gy, flag )

if exist('flag')~=1 flag=0; end;

[sx sy] = size( dummy );
[nx ny] = size( gx );
[xeff, yeff] = bits(sx,sy);

if nx~=xeff | ny~=yeff
  g=zeros(xeff,yeff);
  g(1:nx,1:ny) = gx;  gx=g;
  g(1:nx,1:ny) = gy;  gy=g;
end;

[Gx, Gy] = gradient_essential(dummy, gx, gy);

signal = propagation(Gx,Gy);
% signal = shift(signal);
signal = signal(1:sx,1:sy);

if flag
  figure, imagesc(signal), axis image, colormap gray,
  title('Reconstructed Image'), drawnow;
end;
