%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [P,MM, orMM,  Mx,My,ax,ay,gx,gy] = ...
    exp_sourcediffu(flag, img, entropy, MSM)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% See TEST_POTENTIAL & LAUNCH_PMC
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if exist('flag') ~= 1 flag=0; end;

[sx, sy] = size(img); 
[xeff, yeff] = bits(sx,sy);

flag =1;

fprintf('\nCalcul chromatiquement reduit');
% Calcul chromatiquement reduit
[ax, ay] = derive_MSM_unitary( MSM );
dummy = propagation(ax, ay);
ax = zeros(xeff,yeff);
ax(1:sx,1:sy) = dummy(1:sx,1:sy);


fprintf('\nMethode entropie/chromatiquement reduit');
% methode entropie/chromatiquement reduit
[Gx_entropy, Gy_entropy] = derive_spectral(entropy);
Gx_entropy = Gx_entropy(1:sx,1:sy);
Gy_entropy = Gy_entropy(1:sx,1:sy);
[P,MM, orMM, Mx,My] = ...
    sourcediffu( sx, sy, ax, Gx_entropy, Gy_entropy, flag );

stepsize = .24;
nosteps = 100; 
type='gauss';
k = 25;
g=img;
verb=[1,3,0];
for i=1:nosteps
  
  fprintf('\nDiffusion iter=%d',i);
  fprintf('\nMethode diffuse/chromatiquement reduit');
  % methode diffuse/chromatiquement reduit
  g = pmcstep( img, type, k, stepsize, verb );
  [Gx, Gy] = derive_spectral(g);
  Gx = Gx(1:sx,1:sy);
  Gy = Gy(1:sx,1:sy);
  [P,MM, orMM, Mx,My] = ...
      sourcediffu( sx, sy, ax, Gx, Gy, flag );
  
  fprintf('\nMethode entropie/diffuse');
  % methode entropie/diffuse
  ex = zeros(xeff,yeff);
  ex(1:sx,1:sy) = g;
  [P,MM, orMM, Mx,My] = ...
      sourcediffu( sx, sy, ex, Gx_entropy, Gy_entropy, flag );
  
end;


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [P,MM, orMM, Mx,My] = sourcediffu( sx, sy, ax, gx, gy, flag )
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% See calcula_potential_MSM
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C0=255; CP=0; CM=127;
% C0 = grey level associated to 0
% CM = grey level associated to -1
% CP = grey level associated to 1

EXP_MU=1.;
norma=0.; expon =-EXP_MU;

% See below
% [sx sy] = size( MSM );
[xeff, yeff] = bits(sx,sy);

% 1) Reduced unitary essential vectorial field
% Computes the reduced_from_MSM unitary essential vectorial field
% [ax, ay] = derive_MSM_unitary( MSM );
% % ax=ax(1:sx,1:sy);
% % ay=ay(1:sx,1:sy);
% % 2) 
% dummy = propagation(ax, ay);
    
% figure, imagesc(dummy);

% ax = zeros(xeff,yeff);
% ax(1:sx,1:sy) = dummy(1:sx,1:sy);
% see below
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
