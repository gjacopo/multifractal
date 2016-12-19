%% DERIVE_MSM_UNITARY - Compute the reduced unitary essential.
%
%% Description
% Compute the reduced unitary essential vectorial field which:
%   - is unitary on the MSM, null otherwise,
%   - is perpendicular to the MSM on each pixel of the MSM,
%   - has identical orientation as the original gradient.
%
%% Syntax
%     [Gx, Gy] = derive_msm_unitary( MSM )
%
%% Note 
% The orientation of the gradient is given in the orientated msm.

%% Function implementation
function [Gx, Gy] = derive_msm_unitary( MSM )

% Calcul du gradient unitaire reduit
[ax, ay] = derive_msm( MSM );

% DEBUG
% figure, subplot(1,2,1), imshow(ax,[]);
%         subplot(1,2,2), imshow(ay,[]);
% ENDDEBUG

modo=1;
% remarque: en modo=1, la MSM n'est pas modifiee en fait
% donc pas besoin de la preciser comme argument de sortie
[Gx, Gy] = mask_gradient_msm( msm, ax, ay, ax, ay, modo );
% remarque: les gradients passes ici a la fonction
% mask_gradient_msm sont les memes (ax,ay), donc le produit
% scalaire sigproy (voir mask_gradient_msm.m) sera toujours 
% positif ou nul.

% DEBUG
% figure, subplot(2,2,1), imshow(Gx,[]);
%         subplot(2,2,2), imshow(Gy,[]);
% ENDDEBUG

% Reconstruction de l'image a partir du gradient oriente
Gx = propagation(Gx,Gy);

% DEBUG
% figure, subplot(1,3,1), imshow(Gx,[]);
% return;
% ENDDEBUG

% Derivation du signal reconstruit: the true derive_spectraltive
[Gx,Gy] = derive_spectral( Gx );
% Application du masque de la MSM: mise a 0 en dehors de la MSM
[Gx,Gy] = mask_gradient( MSM, Gx, Gy );
% Le gradient est le gradient de l'image chromatiquement reduite

% DEBUG
%  subplot(2,2,3), imshow(Gx,[]);
%          subplot(2,2,4), imshow(Gy,[]);
% ENDDEBUG
