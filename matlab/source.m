%% SOURCE - Compute the source of an image. 
%
%% Description
% Compute the source of an image, starting from the estimation of the MSM and 
% of the derivatives of the reconstructed image (computed by REDUCED_MSM).
%
%% Syntax 
%      [source, source_vec, Gx, Gy] = source( dummy, gx, gy, silog, flag )
% 
%% Inputs
%   - dummy :chromatically reduced_from_msm image,
%   - gx, gy: gradient vectors of the reconstructed image
%   (remark: these parameters have to be computed by one of the
%   functions reduced_msm or unitary),
%   - silog : if 1, logarithmic representation of source,
%   - flag : if 1, display results.
%
%% Outputs
%   - source: image of source,
%   - source_vec : vectorial representation of source.

%% Function implementation
function [source, source_vec, Gx, Gy] = source( dummy, gx, gy, silog, flag )

if (exist('flag') ~= 1) flag =0; end;

EXP_MU=1.;
norma=0.; expon =-EXP_MU;

% [sx sy] = size( dummy )
[sx sy] = size( gx );
[xeff, yeff] = bits(sx,sy);

ax = zeros(xeff,yeff);
Gx = zeros(xeff,yeff);
Gy = zeros(xeff,yeff);
ax(1:sx,1:sy) = dummy;
Gx(1:sx,1:sy) = gx;
Gy(1:sx,1:sy) = gy;

% filter_spectral: 
%     IFFT( FFT(ax) * f^(-1) ) 
% Multiplication de l'image chromatiquement reduite par f^(-1) dans 
% l'espace des frequences, ou f designe la norme du vecteur frequence.
func = filter_spectral( ax, norma, expon ); % expon=-1, norma=0
% Mise a 0 de la DC-component, de sorte que ax est de moyenne nulle.
% Calcul du gradient unitaire associe a l'image chromatiquement reduite
[ax, ay] = derive_spectral( func );
% Les etapes ci-dessus sont exactement les etapes preliminaires d'estimation 
% du gradient essentiel a partir de l'image chromatiquement reduite (voir
% gradient_essential)

% DEBUG
% figure, imagesc(func(1:sx,1:sy)), colormap gray;
% ENDDEBUG

% Multiplication des derivees par f^(-1) dans l'espace des 
% frequences et mise a 0 de la DC-composante.
Gx = filter_spectral( Gx, norma, expon );
Gy = filter_spectral( Gy, norma, expon ); % expon=-1, norma=0

% DEBUG
% figure, subplot(2,2,1), imshow(Gx(1:sx,1:sy),[]);
%         subplot(2,2,2), imshow(Gy(1:sx,1:sy),[]);
% 	subplot(2,2,3), imshow(ax(1:sx,1:sy),[]);
%         subplot(2,2,4), imshow(ay(1:sx,1:sy),[]);
%	return;
% ENDDEBUG

% Division de (Gx,Gy) par (ax,ay) vus comme les composantes
% de vecteurs complexes. 
% (ax,ay) : gradient unitaire
% (Gx,Gy) : gradient du signal reconstruit
[Gx, Gy] = vec_divide( ax, ay, Gx, Gy );
Gx = Gx(1:sx,1:sy); Gy = Gy(1:sx,1:sy);

% DEBUG
% figure, subplot(2,2,1), imshow(log(Gx),[]);
%         subplot(2,2,2), imshow(log(Gy),[]);
%         subplot(2,2,3), imshow((abs(Gx)<=1e-17)|(abs(Gy)<=1e-17),[]), 
%         title('Scalar source');
% 	Mod = Gx.^2+Gy.^2;
% 	subplot(2,2,3), imshow((Mod<=1e-17),[]);
%         title('Zeros source');
%         Mod = ax.^2 + ay.^2; Mod=Mod(1:sx,1:sy);
%         subplot(2,2,4), imshow((Mod<=1e-17),[]),
% 	title('Poles source');
% return;
% ENDDEBUG

meanx = sum(sum(Gx)) / (sx*sy);
meany = sum(sum(Gy)) / (sx*sy);
fprintf('Source average: (%f,%f):\n', meanx, meany);

normax = std(Gx(:)); normay = std(Gy(:));
fprintf('Source dispersion: %f\n', sqrt(normax.^2+normay.^2));

mmg = min(min(Gx)); mMg = max(max(Gx));
fprintf('Component x: Min.: %f  Max. %f\n',mmg, mMg);
mmg = min(min(Gy)); mMg = max(max(Gy));
fprintf('Component y: Min.: %f  Max. %f\n',mmg, mMg);

% [source, source_vec] = represent_source( Gx, Gy, silog );
source = source_norm( Gx, Gy, silog );
source_vec = source_angle(Gx, Gy);

if flag
  figure, imagesc(source), axis image, colormap gray,
  title('Image of source (log(|\sigma/\sigma_R|))'), drawnow;
  set(gca,'XTickLabel',{''}); set(gca,'YTickLabel',{''});
  figure, imagesc(source_vec), axis image, colormap jet, colorbar,
  title('Orientation of source (\theta = angle(\sigma/\sigma_R))'), drawnow;
  set(gca,'XTickLabel',{''}); set(gca,'YTickLabel',{''});
end;  

