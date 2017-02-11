%% PROPAGATION_CROSS - Reconstruction of a local patch with the micro-wavelet.
%
%%  Syntax
%       RI = propagation_cross( RIx, RIy, Fourier)
%
%% See also
% Related:    
% convolution_cross
% init_cross

%% Function implementation
function RI = propagation_cross( RIx, RIy, Fourier)

j=sqrt(-1);
% Definition of the differential element	
% Naif element:	diff=2.*pi/3.;
diff=sqrt(3.); % id est: 2.*sin(PI/3)

% Cross Fourier Transform
RIx = fourier_cross( RIx, Fourier, -1 ); 
RIy = fourier_cross( RIy, Fourier, -1 ); 

RI = zeros(1,5);
RI(1) = 0.;
RI(2) = (imag(RIx(2)) - j*real(RIx(2))) / diff;
RI(3) = (-imag(RIx(3)) + j*real(RIx(3))) / diff;
RI(4) = (imag(RIy(4)) - j*real(RIy(4))) / diff;
RI(5) = (-imag(RIy(5)) + j*real(RIy(5))) / diff;

% Inverse Cross Fourier Transform
RI = fourier_cross( RI, Fourier, 1 ); 

