%% CONVOLUTION_WAVE - Wavelet convolution.
%
%% Syntax
%     Tmu = convolution_wave( mu, orden, scale, degree, wavetype, flag ) 

%% Function implementation
function Tmu = convolution_wave( mu, orden, scale, degree, ...
				  wavetype, flag )

if (exist('wavetype') ~= 1) wavetype = 'lor'; end;
if (exist('flag') ~= 1) flag =0; end;

[sx sy]=size(mu);
[xeff,yeff] = bits(sx,sy);

fac = sqrt(xeff)*sqrt(yeff);

switch wavetype
   case 'lor'   
    fhandle = @wlorentz;;
    case 'gauss'  
     fhandle = @wgauss; 
%    case ''    
%     fhandle = @; 
end;

% using this psi enables to avoid using fftshift
psi = feval( fhandle, xeff, yeff, orden, scale, degree, flag );
%psi = fspecial('gaussian',[xeff yeff],scale);

FFTpsi = fft2(psi) / fac;
FFTmu = fft2(mu) / fac;

%Tmu = log(abs(ifft2(FFTmu .* FFTpsi) * fac));

Tmu = ifft2(FFTmu .* FFTpsi) * fac;

