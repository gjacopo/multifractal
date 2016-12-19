%% VEC_MULTIPLY - Complex dot multiplication.
%
%% Description
% Complex multiplication of two vectorial fields seen as
% components in the complex plane.
%
%% Syntax :
%     [Gx, Gy] = vec_multiply( ax, ay, gx, gy)
%
% Outputs
%     - Gx = (ax.*gx - ay.*gy),
%     - Gy = (ax.*gy + ay.*gx). 
%
%% See also  
% Related:
% <|gradient_essential|>

function [Gx, Gy] = vec_multiply( ax, ay, gx, gy) 

Gx = (ax.*gx - ay.*gy) ;
Gy = (ax.*gy + ay.*gx) ;

