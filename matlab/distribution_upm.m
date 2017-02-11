%% DISTRIBUTION_UPM - Compute the values and the distribution of singularity 
% exponents.
%
%% Syntax
%      Disp = distribution_upm( img, gx, gy )
% 
%% Inputs
%     - img: original image,
%     - gx, gy: derivatives of img (computed with derive_spectral).
%
%% Outputs
%     - Disp: matrice of singularity exponents for exch pixel
%       of the image.
%
%% See also
% Related:    
% upm

%% Function implementation
function Disp = distribution_upm( img, gx, gy )

width=1;

[sx sy] = size(img);
[xeff yeff]= size(gx); % puissances de 2
% sx: nb lignes, sy: nb colonnes
% on duplique les bords de img
img = duplicate(img,1);

% Complex matrix for Cross Fourier Transform: Fourier{1}
% Complex matrix for Inverse Cross Fourier Transform: Fourier{2}
Fourier=init_cross(0);

Errx = zeros(xeff,yeff);
Erry = zeros(xeff,yeff);
Disp = zeros(xeff,yeff);

for iy=1:sy
  for ix=1:sx
    parche = clique( img, ix+1, iy+1 );
    [errx, erry] = convolution_cross(parche, Fourier);
    Errx(ix,iy) = errx;
    Erry(ix,iy) = erry;
  end;
end;

filt = ones(2*width+1,2*width+1); 

Gxy2 = gx.^2 + gy.^2;
Gxy2 = duplicate(Gxy2,width);
Proy = conv2( Gxy2, filt,'same');

% Somme locale des erreurs dans la composantes x
Errx= duplicate(Errx,width);
ErrxS = conv2( Errx, filt,'same');
% Somme locale des erreurs dans la composantes y
Erry= duplicate(Erry,width);
ErryS = conv2( Erry, filt,'same');
Mod = ErrxS.*Errx + ErryS.*Erry;

% Mod=Mod(1+width:sx+width,1+width:sy+width)
% Proy=Proy(1+width:sx+width,1+width:sy+width)

Disp = (Proy~=0.) .* sqrt(abs(Mod./Proy)) + (Proy==0.).*1;
Disp = Disp(1+width:sx+width,1+width:sy+width);


