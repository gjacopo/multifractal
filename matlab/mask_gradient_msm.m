%% MASK_GRADIENT_MSM - Compute an orientated MSM.
%
%% Description
% Compute an orientated MSM (or not, depending on the orientation of the gradient 
% over MSM) and signed (or not) gradient's vectors according to this orientation.
%
%% Syntax 
%     [Gx, Gy, MSM] = mask_gradient_msm( msm, ax, ay, gx, gy, modo )
%
%% Notes
%   1) in reduced_msm: [gx,gy]=derive_spectral(img): gradient de 
%      l'image originale;
%      in derive_msm_unitary: [gx,gy]=derive_msm(msm);
%   2) we will always use the gradient (ax,ay) as the result of:  
%            [ax, ay] = derive_msm( msm ); 
%      but this is done elsewhere to avoid recurrent reprocessing.
%
%% See also
% Related:    
% derive_msm_unitary
% reduced_msm

%% Function implementation
function [Gx, Gy, MSM] = mask_gradient_msm( msm, ax, ay, gx, gy, modo )

C0=255; CP=0; CM=127;
% C0 = grey level associated to 0
% CM = grey level associated to -1
% CP = grey level associated to 1

[sx sy] = size(msm);
[xeff, yeff] = bits(sx,sy);

Gx = zeros(xeff, yeff);
Gy = zeros(xeff, yeff);

% not done here: done upper in the calling functions
% [ax, ay] = derive_msm( msm );

% Calcul du produit scalaire des vecteurs (ax,ay) et (gx,gy)
% sur toute l'image
sigproy = ax.*gx + ay.*gy; % de dim. [xeff,yeff]
sigproy = (msm~=C0) .* sigproy(1:sx,1:sy); % + (msm==C0)*0
% sigproy renseigne sur l'orientation du gradient (gx,gy) par rapport
% au gradient unitaire (ax,ay) sur les pixels de la MSM sur les
% pixels de la MSM (nul en dehors). Suivant le signe de sigproy:
%   - sigproy > 0 : les 2 vecteurs (gx,gy) et (ax,ay) ont la meme
%     orientation en ce point (dans un meme quadrant)
%   - sigproy < 0 : les 2 vecteurs ont des orientations opposees

gxaux = gx(1:sx,1:sy);
gyaux = gy(1:sx,1:sy);

if(modo==0) 
  % same as mask_gradient for Gx and Gy
  Gx(1:sx,1:sy) = (msm~=C0).*gxaux;
  Gy(1:sx,1:sy) = (msm~=C0).*gyaux;
  % give a sign to the MSM according to the orientation of the gradient
  MSM = (sigproy<0) .* CM + (sigproy>=0) .* msm;
else 
  mask =  ((sigproy>=0.)&(msm==CM)) | ((sigproy<0.)&(msm==CP)); % (msm~=C0) .*
  mask = mask + (msm==C0).*(-1);
  % Orientate the gradient
  Gx(1:sx,1:sy) = (mask==1) .* (-gxaux) + (mask==0) .* (gxaux);
  Gy(1:sx,1:sy) = (mask==1) .* (-gyaux) + (mask==0) .* (gyaux);
  % The MSM is unchanged
  MSM = msm;
end;


