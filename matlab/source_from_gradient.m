%% SOURCE_FROM_GRADIENT - Compute the source of a field.
%
%% Description
% Compute the source of a field which is known through its MSM and its 
% gradient values.
%
%% Syntax 
%     [MM, orMM, Mx, My] = source_from_gradient( MSM, gx, gy, flag )
% 
%% Inputs
%     - MSM : orientated MSM,
%     - gx, gy : gradient vector field as if they were outputs
%       of the function reduced_MSM,
%     - flag : if 1, display some intermedary results,
%     - MM : lognorm of the source,
%     - orMM : orientation of the source,
%     - Mx, My : gradient vector field of source.
% 
%% See also
% Related:     
% run_fractal_entropy
% source

%% Function implementation
function [MM, orMM, Mx, My] = source_from_gradient( MSM, gx, gy, flag )

C0=255; CP=0; CM=127;
% C0 = grey level associated to 0
% CM = grey level associated to -1
% CP = grey level associated to 1

EXP_MU=1.;
norma=0.; expon =-EXP_MU;

[sx sy] = size( MSM );
[xeff, yeff] = bits(sx,sy);

% Computes the reduced_from_MSM unitary essential vectorial field
% and then computes the chromatically reduced_from_MSM image.
dummy =  dummy( MSM );    
    
% figure, imagesc(dummy);

ax = zeros(xeff,yeff);
ax(1:sx,1:sy) = dummy(1:sx,1:sy);
% filter_spectral: %     IFFT( FFT(ax * f^(-1)) ) 
% Multiplication de l'image chromatiquement reduite par f^(-1) dans 
% l'espace des frequences, ou f designe la norme du vecteur frequence.
func = filter_spectral( ax, norma, expon ); % expon=-1, norma=0
% Mise a 0 de la DC-component, de sorte que ax est de moyenne nulle.
% Calcul des derivees
[ax, ay] = derive_spectral( func );

Gx=zeros(xeff,yeff); Gy=Gx;
Gx(1:sx,1:sy) = gx;
Gy(1:sx,1:sy) = gy;
gx = filter_spectral( Gx, norma, expon );
gy = filter_spectral( Gy, norma, expon ); % expon=-1, norma=0

if flag
  figure, subplot(2,2,1), imshow(ax(1:sx,1:sy),[]), colormap(gray), ...
      title('ax');
  subplot(2,2,2), imshow(ay(1:sx,1:sy),[]), colormap(gray), ...
      title('ay');
  subplot(2,2,3), imshow(gx(1:sx,1:sy),[]), colormap(gray), ...
      title('gx');
  subplot(2,2,4), imshow(gy(1:sx,1:sy),[]), colormap(gray), ...
      title('gy');
  suptitle('Unitary (top) and reconstructed (bottom) gradients');
end;

[Mx, My] = vec_divide( ax, ay, gx, gy );
Mx = Mx(1:sx,1:sy);
My = My(1:sx,1:sy);
MM = source_norm(Mx, My,1);
orMM = source_angle(Mx, My);

if flag
  figure, imagesc(MM), colormap(gray), axis image, title('Sources');
  set(gca,'XTickLabel',{''}); set(gca,'YTickLabel',{''});
  figure, imagesc(orMM), axis image, colormap jet, colorbar, 
  title('Sources orientation');
  set(gca,'XTickLabel',{''}); set(gca,'YTickLabel',{''});
end;
