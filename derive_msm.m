% DERIVE_MSM - Transform the Hausdorff measure of the MSM in the frequency domain.
%
%
%% Description
% The function computes:
%     IFFT( FFT(\delta_{msm} * f^(-1)) ) 
% where \delta_{msm} stands for the Hausdorff measure restricted
% to the MSM and f stands for the norm of the frequency vector.
% Normalizes the result (unitary vector).
%
%% Syntax :
%       [ax, ay] = derive_msm( MSM )
%
%% Outputs
% Returns: ax, ay of dimension [xeff,yeff]=bits(MSM)        
%
%% See also
% Related:    
% mask_gradient_msm
% derive_msm_unitary

%% Function implementation
function [ax, ay] = derive_msm( msm )

C0=255; CP=0; CM=127;
% C0 = grey level associated to 0
% CM = grey level associated to -1  
% CP = grey level associated to 1

[sx sy] = size(msm);
[xeff, yeff] = bits(sx,sy);

ax = zeros(xeff,yeff);

% pinta_suave
% Calcul du masque de la MSM: ax = \delta_{msm}
ax(1:sx,1:sy) = (msm~=C0);
% Multiplication du masque de la MSM par f^(-1) dans l'espace des
% frequences, ou f designe la norme du vecteur frequence.
ax = filter_spectral( ax, 0., -1. );
% Mise a 0 de la DC-component, de sorte que ax est de moyenne nulle.
% fin pinta_suave

% DEBUG
% sum(sum(ax))
% figure, imagesc(ax),colormap gray
% ENDDEBUG

% Compute the derivatives
[ax, ay] = derive_spectral( ax );
% ... and normalize
[ax, ay] = vec_normalise(ax, ay);

% DEBUG
% figure, subplot(1,2,1), imshow(ax,[]);
%         subplot(1,2,1), imshow(ay,[],'turesize');
% ENDDEBUG
