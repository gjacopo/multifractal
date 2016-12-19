%% CONVOLUTION_CROSS - Convolution with the micro-wavelet.
%
%% Syntax
%         [errx, erry] = convolution_cross(parche, Fourier)
%
%% Inputs
%    - parche : patch (neighborhood) extracted from the signal,
%    - Fourier : micro-wavelet defined in Direct and in Fourier
%      space for decomposition/reconstruction (see init_cross).

%% Function implementation
function [errx, erry] = convolution_cross(parche, Fourier)

diff =sqrt(3);
j = sqrt(-1);

% Somme du signal dans le voisinage considere 
mediap=sum(parche)/3.;

% Ponderation locale du voisinage
parche(1) = parche(1) + mediap;
for ip=2:5
  parche(ip) = parche(ip) - mediap;
end;
% On a alors :  
%     parche[0] = - sum_{ip=1}^4 parche[ip] 
% =>  sum_{ip=0}^4 parche[ip] = 0    

% Calcul du signal derive avec la micro-ondelette
RI = fourier_cross( parche, Fourier, -1 );

% Composante en x
parche(1) = 0.;
parche(2) = -diff * imag(RI(2)) ...
    + j * diff * real(RI(2));
parche(3) = diff*imag(RI(3)) ...
    - j * diff * real(RI(3));

% Composante en y 
RIgy = zeros(1,5);
RIgy(4) = -diff * imag(RI(4)) ...
	 + j * diff * real(RI(4));
RIgy(5) = diff * imag(RI(5))...
	 - j * diff * real(RI(5));

parche(4) = 0.;
parche(5) = 0.;

errx = real(parche(2)) + real(parche(3));
erry = real(RIgy(4)) + real(RIgy(5));

RIx = fourier_cross( parche, Fourier, 1 );
RIy = fourier_cross( RIgy, Fourier, 1 );

RIx(1) = 0.;
RIy(1) = 0.;
parche = fourier_cross( RIx, Fourier, -1 );
RIgy = fourier_cross( RIy, Fourier, -1 );

errx = errx + real(parche(1));
erry = erry + real(RIgy(1));
