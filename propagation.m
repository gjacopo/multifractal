%% PROPAGATION - Reconstruct a signal from the essential density and the 
% universal propagator.
%
%% Description
% Reconstruct a signal c from the essential (vector-valued) density and the 
% universal propagator g:
%         c = g ** v
% where ** stands for the convolution product (being understood as a scalar 
% products of vectors).
%
%% Syntax
%         rec = propagation(Gx,Gy)
%
%% Note
% The propagator g is given in the Fourier space by:
%          F(g) = i*f / |f|^2
% where i is the imaginary unit (i=sqrt(-1)) and f stands for 
% the frequency vector. So to retrieve the signal, computations 
% are made in the Fourier space.
%
%% See also
% Related:    
% derive_msm_unitary 
% upm
% gradient_essential
% reduced_unitary

%% Function implementation
function rec = propagation(Gx,Gy)

j=sqrt(-1);

[xeff yeff] = size(Gx);
fac = sqrt(xeff)*sqrt(yeff);

% Direct Fourier Transform of the signal in Gx
FFTgx = fft2(Gx) / fac;
% ce qui equivaut a: FFTgx = fft(fft(gx/sqrt(sx)).'/sqrt(sy)).';
% Direct Fourier Transform of the signal in Gy
FFTgy = fft2(Gy) / fac;

[X Y] = meshgrid( [0:yeff/2-1,0,-yeff/2+1:-1]/yeff, ...
		  [0:xeff/2-1,0,-xeff/2+1:-1]/xeff);
% Vecteur frequence
dX = sin(pi.*X);
dY = sin(pi.*Y);
% or: dX =X; dY =Y;

prefm = dX.^2 + dY.^2;
% probleme de la division par zero quand dX=dY=0
prefm = (prefm>1e-30) .* prefm + (prefm<=1e-30).*1;

% The scalar product with the frequency vector f is computed
% and, simultaneously, the complex number so obtained is multiplied by
% the imaginary unit i=sqrt(-1), wich correspond to a rotation:
%     A = a_0 + a_1*i   =>   B = A*i = -a_1 +a_0*i         
gxR = (dX.*imag(FFTgx) + dY.*imag(FFTgy)) ./ prefm;
gxI = - (dX.*real(FFTgx) + dY.*real(FFTgy)) ./ prefm;

% Inverse Fourier transform 
rec = real(ifft2( gxR+j*gxI )) * fac;
