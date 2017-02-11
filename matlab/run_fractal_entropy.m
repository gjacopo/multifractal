%% RUN_FRACTAL_ENTROPY - Compute the source associated with the 
% gradient field of entropy. 
%
%% Syntax
%   [MM, orMM, Mx, My, entropy, MSM] = ...
%         run_fractal_entropy( [ img, entropy, MSM, flag] )

%% Function implementation
function [MM, orMM, Mx, My, entropy, MSM] = ...
    run_fractal_entropy( varargin )

% if exist('flag') ~= 1 flag=0; end;

%% Load or compute the needed images
if nargin <= 1 | (nargin ==2 & size(varargin{2},1)==1 & size(varargin{2},2)==1)
  img = varargin{1};
  % Computes the multiscale entropy of the image.
  if flag  fprintf('\nComputing the entropy'); end;
  mode=-2; width=10;
  entropy = margEntropy2D( img, mode, width );
  fprintf('\n Computing the MSM and the chromatically reduced_from_msm image');
  % Computes the MSM and the chromatically reduced_from_msm image of the
  % image (associated with the true gradient of image).
  upm_dens=0.45; upm_thres=1.;
  [MSM, dummy, Gx, Gy] = ...
      reduced_msm(img, upm_dens, upm_thres, flag);
  if (nargin ==2) flag=varargin{2}; end;
  
elseif nargin <= 2 | (nargin ==3 & size(varargin{3},1)==1 & size(varargin{3},2)==1)
  img = varargin{1};
  entropy = varargin{2};
  % Computes the MSM and the chromatically reduced_from_msm image of the
  % image (associated with the true gradient of image).
  fprintf('\n Computing the MSM and the chromatically reduced_from_msm image');
  upm_dens=0.45; upm_thres=1.;
  [MSM, dummy, Gx, Gy] = ...
      reduced_msm(img, upm_dens, upm_thres, flag);
  if (nargin ==3) flag=varargin{3}; end;  

elseif nargin <=4
  img = varargin{1};
  entropy = varargin{2};
  MSM = varargin{3};
  if (nargin ==4) flag=varargin{4}; else flag =0;  end;
end;

if flag
  figure, subplot(1,2,1), imshow(img,[]), title('Original image');
  subplot(1,2,2), imshow(entropy,[]), title('Multiscale entropy');
  drawnow;    
end;

%% Analyze the mean value of entropy over the MSM
[sx, sy] = size(img);
plot_on_msm( MSM, entropy, flag);

%% Defines the gradient vector field associated with entropy 
% Computes the gradient of entropy
[Gx, Gy] = derive_spectral(entropy);
Gx=Gx(1:sx,1:sy);
Gy=Gy(1:sx,1:sy);
% Computes the orientation of the entropy field
output = angle(Gx, Gy);

if flag 
  figure, subplot(1,4,1), imshow(Gx,[]);
  subplot(1,4,2), imshow(Gy,[]);
  subplot(1,4,3), imshow(log(Gx.^2+Gy.^2), []);
  suptitle('Gradients of entropy and lognorm');
  figure, imagesc(output), colormap jet, colorbar, 
  title ('Orientation of entropy field');
end;	

%% Main computation
% Computes the source of the field defined by the gradient of
% entropy over the MSM of the image
[MM, orMM, Mx,My] = ...
    source_from_gradient( MSM, Gx, Gy , flag );
