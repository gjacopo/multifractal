%% WLORENTZ - Compute a lorentzian wavelet at a given scale.
%
%% Syntax
%     psi = wlorentz( xeff, yeff, orden, scale, degree [, flag] )
%
%% See also 
% Related: 
% <|convolution_wave|>

%% Function implementation
function psi = wlorentz( xeff, yeff, orden, scale, degree, flag )

if (exist('flag','var') ~= 1), flag =0; end;

[X Y] = meshgrid( [0:xeff/2-1,-xeff/2:-1], ...
		  [0:yeff/2-1,-yeff/2:-1] );

a0=scale;
R = sqrt(X.*X + Y.*Y) / scale;

if degree==0
  pref = ones(xeff,yeff);
elseif degree==1
    pref = -R;
elseif degree==2
    pref = (2.*orden+1.) .* R.^2 -1.;
elseif degree==3
    pref = -R .* ((2.*orden+1.) .* R.^2 -3.);
elseif degree==4
    pref = (2.*orden+1.)*(2.*orden+3.) .* R.^4 ...
	   -6.*(2.*orden+3.) .* R.^2 +3.;
else
  pref = ones(xeff,yeff);  
end;

psi = pref .* (1+R.^2) .^ (-(orden+degree)) ./ (a0*a0);

if flag
  [X Y] = meshgrid( 1:xeff, 1:yeff );
  figure, mesh(X,Y,psi);
end;
