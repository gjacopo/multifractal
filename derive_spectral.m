%% DERIVE_SPECTRAL - Compute image derivatives of an image in frequency domain 
% using the micro-wavelet.
%
%% Syntax
%       [ax, ay] = derive_spectral( f )
%
%% Note
% The derive_spectraltive is omputed through the Fourier formula: if FT stands 
% for the Fourier transform and ' for the derivation, we have for any integrable 
% function f:
%            FT{f'}(y) = iy FT{f}(y)
% where y is the frequency and i the imaginary unit, so the derivative of f can 
% be expressed wrt its Fourier transform FT{f} by:
%            f'(x) = IFT{ y->iy FT{f}(y)}
% where IFT stands for the inverse Fourier transform.
%
%% See also
% Related: 
% derive_msm   
% derive_msm_unitary 
% msm
% gradient_essential
% source

%% Function implementation
function [ax, ay] = derive_spectral( f )

j=sqrt(-1);

[sx sy] = size( f ); % sy: nbre de lignes, sx: nbre de colonnes
[xeff yeff] = bits(sx,sy);
% remarque: xeff et yeff seront toujours des puissances de 2
fac = sqrt(xeff)*sqrt(yeff);

g=zeros(xeff, yeff);
g(1:sx,1:sy) = f;
% Direct Fourier Transform of the signal in gx
FFTgx = fft2(g) / fac;
% ce qui equivaut a: FFTgx = fft(fft(gx/sqrt(xeff)).'/sqrt(yeff)).';

[X Y] = meshgrid( [0:yeff/2-1,0,-yeff/2+1:-1]/yeff, ...
		  [0:xeff/2-1,0,-xeff/2+1:-1]/xeff);
% ou de maniere equivalente:
% [X Y] = meshgrid([0:yeff/2,-yeff/2+1:-1]/yeff,[0:xeff/2,-xeff/2+1:-1]/sx)
% X(floor(xeff/2)+1,:) = 0.;
% Y(:,floor(yeff/2)+1) = 0.;

% Vecteur frequence
dX = sin(pi.*X);
dY = sin(pi.*Y);
% or: dX =X; dY =Y;

% Inverse Fourier Transform
ax = real(ifft2( - dX.*imag(FFTgx) +  j*dX.*real(FFTgx) )) * fac;
ay = real(ifft2( - dY.*imag(FFTgx) +  j*dY.*real(FFTgx) )) * fac;
