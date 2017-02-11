%% WGAUSS - Compute a lorentzian wavelet at a given scale.
%
%% Syntax
%     [psi,sc] = wgauss( xeff, yeff, scale [, flag] )
%
%% See also
% Related:
% <|convolution_wave|>

%% Function implementation
function [psi,sc] = wgauss( xeff, yeff, orden, scale, degree, flag )   %#ok

if (exist('flag','var') ~= 1),  flag =0; end;

[X Y] = meshgrid( [0:xeff/2-1,-xeff/2:-1], ...
		  [0:yeff/2-1,-yeff/2:-1] );
%[X Y] = meshgrid( [-xeff/2:-1, 0:xeff/2-1], ...
%    [-yeff/2:-1,0:yeff/2-1]  );
% the previous meshgrid avoids the futher call to fftshift
% in convolution processes...

R = (X.*X + Y.*Y);

psi = exp(- R / (2*scale^2)); % / scale^2;

s=sum(abs(psi(:)))
if s>1.e-17
    psi = psi/s;
end
norma = s/(xeff*yeff);
sc = sqrt(norma/pi);

sum(psi(:))
[min(psi(:)) max(psi(:))]

if flag
  [X Y] = meshgrid( 1:xeff, 1:yeff );
  figure, mesh(X,Y,psi);
end;
