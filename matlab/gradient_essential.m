%% GRADIENT_ESSENTIAL - Compute the essential gradient from the chromatically 
% reduced image.
%
%% Description
% Compute the essential gradient field from the chromatically reduced image. 
% The essential gradient field is the gradient restricted to the MSM.
%
%% Syntax
%     [Gx, Gy] = gradient_essential(dummy, gx, gy)
%
%% Inputs
%   - dummy :chromatically reduced_from_msm image (computed by 
%     reduced_msm or unitary),
%   - gx, gy: gradient images computed by source.
%
%% See also
% Related:    
% reconstruction

%% Function implementation
function [Gx, Gy] = gradient_essential(dummy, gx, gy)

EXP_MU=1.;
norma=0.; expon=EXP_MU;

[sx sy] = size( dummy );
[xeff, yeff] = bits(sx,sy);

% done upper, in reconstruction
% [nx ny] = size( gx );
% if nx~=xeff | ny~=yeff
%   g=zeros(xeff,yeff);
%   g(1:nx,1:ny) = gx;  gx=g;
%   g(1:nx,1:ny) = gy;  gy=g;
% end;

ax = zeros(xeff,yeff);
ax(1:sx,1:sy) = dummy;
func = filter_spectral( ax, norma, -expon );
[ax, ay] = derive_spectral( func );
% Les operations ci-dessous sont deja realisees dans le processus
% de construction des source de luminance (voir source)

[Gx, Gy] = vec_multiply(ax, ay, gx, gy);
Gx = propagation(Gx, Gy);

func = filter_spectral( Gx, norma, expon );
[Gx, Gy] = derive_spectral( func );
