%% ANGLE - Compute the orientation of any vector.
%
%% Syntax     
%           output = angle(gx, gy )
%
%% Note
% The output lies in the four quadrant arctangent of the real parts of the
% elements of gx and gy:        0 <= output <= 2*pi

%% Function implementation
function output = angle(gx, gy )

% ou : 
output=atan2(gy,gx);
% , alors -pi <= output <= pi
output = (output +2*pi).*(output<0)  + ...
	 output .* (output>=0);

%  P = atan2(Y,X) returns an array P the same size as X and Y 
% containing the element-by-element, four-quadrant inverse tangent
% (arctangent) of the real parts of Y and X. 
% Any imaginary parts are ignored.
% Elements of P lie in the closed interval [-pi,pi]. The specific
% quadrant is determined by sign(Y) and sign(X). 
% This contrasts with the result of atan(Y/X), which is limited to
% the interval[-pi/2, pi/2]


return;

m = (abs(gx)>1e-30);
gx = m .* gx + ~m;
output = atan(gy ./ gx) + pi * (gx<0);

output = output + 2*pi * (output<0);
output = output -2*pi * (output>2*pi);

output = m .* output ...
	 + ((m)&(gy>0)) .* pi/2. + ((~m)&(~(gy>0))) .* 3.*pi/2.;

