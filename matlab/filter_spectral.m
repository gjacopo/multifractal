%% FILTER_SPECTRAL - Perform filtering in the frequency domain.
%
%% Description
% Transformation in the frequency domain: scalar product with the norm of the
% frequency vector in the frequency domain.
%
%% Syntax
%      func = filter_spectral( funcion, norma, expon )
% 
%% Inputs
%     - funcion: matrice to be transformed in the frequency 
%       domain,
%     - norma: coefficient multiplication of the DC-component 
%       (0-frequency),
%     - expon: exponent of the power of f wich will be used
%       in the multiplication of the FT of funcion, where f
%       stands for the frequency vector.
%
% Notes
% The function computes:
%     IFT( FT(func) * f^expon) )
% namely:
%       / FT(funcion) * f^expon   for non-nul frequencies,
%       \ FT(funcion) * norma     for DC-component,
% where f=(2sin(pi.x), 2sin(pi.y)) stands for the norm of the 
% frequency vector, , and eventually puts the DC-component to 0,
% and then computes the IFT.
%
%% See also
% Related:    
% derive_msm
% source
% gradient_essential

%% Function implementation
function func = filter_spectral( funcion, norma, expon )

% funcion doit etre de dimensions (puissance de 2)

[yeff xeff] = size(funcion);
fac = sqrt(xeff)*sqrt(yeff);

% Direct Fourier Transform
FFTfunc = fft2(funcion) / fac;

[X Y] = meshgrid( [0:xeff/2,-xeff/2+1:-1]/xeff, ...
		  [0:yeff/2,-yeff/2+1:-1]/yeff );
dX = 2.*sin(pi.*X);
dY = 2.*sin(pi.*Y);
% or: dX =X; dY =Y;

f = sqrt(dX.^2 + dY.^2);
f = (f>1e-30) .* f.^expon;

aux = FFTfunc(1,1);
FFTfunc = FFTfunc .* f;
FFTfunc(1,1) = aux * norma; 

func = real( ifft2(FFTfunc) ) * fac;
