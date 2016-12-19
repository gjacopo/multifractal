%% INIT_CROSS - Initialise the micro-wavelet for analysis/reconstruction. 
%
%% Syntax
%      Fou = init_cross([flag])
%
%% Input
% Parameter: if flag==1, display the microwavelet
%
%% Note
% Matrix for Cross Fourier Transform (real and imaginary parts)	
%
%---------------------------------|--------------------------|
%            Fou[1] (real)        |     Fou[1] (imag.)       | 
% --------------------------------|--------------------------|
%                                 |                          |
%   1/3  1/3   1/3   1/3   1/3    |    a  0   0   0   0      |
%   1/3 -1/6  -1/6   1/3   1/3    |    0  a  -a   0   0      |
%   1/3 -1/6  -1/6   1/3   1/3    |    0 -a   a   0   0      |
%   1/3  1/3   1/3  -1/6  -1/6    |    0  0   0   a  -a      |
%   1/3  1/3   1/3  -1/6  -1/6    |    0  0   0  -a   a      |  
%                                 |                          |
%                                 | a=sqrt(3.)/6.            |
% --------------------------------|--------------------------|
% 
% Inversa matrix for Cross Fourier Transform (real and imaginary parts)
%---------------------------------|--------------------------|
%          Fou[3] (real)          |    Fou[3] (imag.)        |    
%---------------------------------|--------------------------|
%                                 |                          |
%  -1   1     1     1     1       |    0  0   0   0   0      |
%   1  -1/2  -1/2   0     0       |    0 -b   b   0   0      |
%   1  -1/2  -1/2   0     0       |    0  b  -b   0   0      | 
%   1   0     0    -1/2  -1/2     |    0  0   0  -b   b      |
%   1   0     0    -1/2  -1/2     |    0  0   0   b  -b      |
%                                 |                          |
%                                 | b=sqrt(3.)/2.            |
% --------------------------------|--------------------------|

%% Function implementation
function Fou = init_cross(flag)

if (exist('flag') ~= 1) flag=0; end;
  
Fou = cell(2,5,5);
j = sqrt(-1);

Fou{1} = ones(5,5) * 1/3.;
%Fou{2} = zeros(5,5);
Fou{2} = zeros(5,5);
%Fou{4} = zeros(5,5);

for ipx=2:5
  Fou{1}(ipx,ipx) = -1/6. + j * sqrt(3.)/6.;
  % Fou{2}(ipx,ipx) = sqrt(3.)/6.;
  Fou{2}(ipx,ipx) = -1/2. - j * sqrt(3.)/2.;
  Fou{2}(ipx,1) = 1.;
  Fou{2}(1,ipx) = 1.;
  % Fou{4}(ipx,ipx) = -sqrt(3.)/2.;
end;


Fou{1}(2,3) = -1/6. - j * sqrt(3.)/6.;
Fou{1}(3,2) = -1/6. - j * sqrt(3.)/6.;
Fou{1}(4,5) = -1/6. - j * sqrt(3.)/6.;
Fou{1}(5,4) = -1/6. - j * sqrt(3.)/6.;

% Fou{2}(2,3) = -sqrt(3.)/6.;
% Fou{2}(3,2) = -sqrt(3.)/6.;
% Fou{2}(4,5) = -sqrt(3.)/6.;
% Fou{2}(5,4) = -sqrt(3.)/6.;

Fou{2}(1,1) = -1;


Fou{2}(2,3) = -1/2. + j * sqrt(3.)/2.;
Fou{2}(3,2) = -1/2. + j * sqrt(3.)/2.;
Fou{2}(4,5) = -1/2. + j * sqrt(3.)/2.;
Fou{2}(5,4) = -1/2. + j * sqrt(3.)/2.;

% Fou{4}(2,3) = sqrt(3.)/2.;
% Fou{4}(3,2) = sqrt(3.)/2.;
% Fou{4}(4,5) = sqrt(3.)/2.;
% Fou{4}(5,4) = sqrt(3.)/2.;


if flag 
  [X, Y] = meshgrid(1:5);
  for i=1:2
    figure, surf( X, Y, real(Fou{i}), 'FaceColor','interp', 'EdgeColor','none',...
		  'FaceLighting','phong'); colormap gray;
    axis tight, camlight left;
  
    figure, imagesc(real(Fou{i})), axis image, colormap gray
    set(gca,'XTickLabel',{'','1','','2','','3','','4','','5'});
    set(gca,'YTickLabel',{'','1','','2','','3','','4','','5'});

    
    figure, surf( X, Y, imag(Fou{i}), 'FaceColor','interp', 'EdgeColor','none',...
		  'FaceLighting','phong'); colormap gray;
    axis tight, camlight left;
  
    figure, imagesc(imag(Fou{i})), axis image, colormap gray
    set(gca,'XTickLabel',{'','1','','2','','3','','4','','5'});
    set(gca,'YTickLabel',{'','1','','2','','3','','4','','5'});
  
  end;
end;
