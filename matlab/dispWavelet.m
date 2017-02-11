%% DISPWAVELET - Displays 2D Lorentzian and Gaussian wavelets.
%
%% Syntax   
%         dispWavelet( dim, type, gamma )
%
%% Inputs
%       - dim : 1 or 2, dimension of the space,
%       - type : 'lorentz' or 'gauss',
%       - gamma : coefficient of gaussian or lorentzian wavelet

%% Function implementation
function dispWavelet( dim, type,  gamma )

if exist('gamma') ~= 1 | gamma==0
  gamma=1; 
end;

width = 50; % pas de discretisation
winsize=2*width+1;

%% Graphe en 1D
if dim == 1
x=(-width:width)/winsize;
r = x.^2;

if strcmp(type,'lorentz')
  psi = 1. ./ ((1+r).^gamma);
elseif strcmp(type,'gauss')
  psi = exp(-r/(2.*gamma^2))./gamma^2;  
end;
figure, plot(x,psi);

%% Graphe en 2D
elseif dim == 2
  % Calcul de |x|^\alpha
[x,y]=meshgrid(-width:width,-width:width);
x=x/winsize;
y=y/winsize;
r = x.^2 + y.^2;
r(width+1,width+1);

if strcmp(type,'lorentz')
  psi = 1. ./ ((1+r).^gamma);
elseif strcmp(type,'gauss')
  psi = exp(-r./(2.*gamma^2))./gamma^2;  

end;

% figure, meshc(x,y,psi), axis equal, axis square, shading interp
figure, surf( x, y, psi,'FaceColor','interp', 'EdgeColor','none',...
    'FaceLighting','phong');  shading interp
axis square,axis tight, camlight left;

end;
