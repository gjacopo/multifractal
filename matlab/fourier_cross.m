%% FOURIER_CROSS - Cross Fourier Transform or Inverse Cross Fourier Transform 
% with the micro-wavelet.
%
%% Syntax
%        RI = fourier_cross( RIf, Fourier, sign ) 
%
%% See also
% Related:    
% propagation_cross
% convolution_cross

%% Function implementation
function RI = fourier_cross( RIf, Fourier, sign ) 

if sign<0 iFou = 1; else iFou=2; end;

RI = zeros(1,5);

% Multiplication complexe pour la convolution de Rf avec
% Fourier{ifou} 
for ipx=1:5
  aux = Fourier{iFou}(ipx,:) .* RIf;
  RI(ipx) = sum(aux);
end;


