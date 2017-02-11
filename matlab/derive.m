%% DERIVE - Compute naive derivatives of an image. 
%
%% Syntax
%       [ay, ax] = derive( img )

%% Function implementation
function [ay, ax] = derive( img )

[sx sy] = size(img); 
dup = duplicate(img, 1);
% ou: dup(2:sx+1,2:sy+1) = img

ax=dup(3:sx+2,2:sy+1)-dup(2:sx+1,2:sy+1);
ay=dup(2:sx+1,3:sy+2)-dup(2:sx+1,2:sy+1);
