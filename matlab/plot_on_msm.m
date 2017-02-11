%% PLOT_ON_MSM - Compute the distribution of any feature on the MSM. 
%
%% Syntax
%    [f_on_msm, f_out_msm] = plot_on_msm( msm, feature, flag )
%
%% Note
% The MSM is used as a mask. Values in the mask of MSM must be C0 (not in the 
% MSM), CP or CM (in the MSM).

%% Function implementation
function [f_on_msm, f_out_msm] = plot_on_msm( msm, feature, flag )

if exist('flag') ~=1 flag=0; end;

C0=255.;
% C0 = grey level associated to 0: pixels not belonging to the MSM

[sx sy]=size(msm);
[nx ny]=size(feature);
if nx~=sx | ny~=sy
  fprintf('\n!in plot_on_msm: input matrices do not have identical dimensions!');
  sx=min(sx,nx); sy=min(ny,sy);
  msm=msm(1:sx,1:sy); feature=feature(1:sx,1:sy);
end;  
  
Nhisto=1000;

I=find(msm~=C0);
dens=size(I) / (sx*sy);
f_on_msm = feature(I);
if flag fprintf('\nMSM at density %f',dens); end;
h_on_msm = hist(f_on_msm,Nhisto);
% remark: feature(I) is a vector
norerror = std(f_on_msm);
fprintf('\nDispersion of the feature over the MSM: %f',norerror);

I=find(msm==C0);
f_out_msm = feature(I);
h_out_msm = hist(f_out_msm,Nhisto);
norerror = std(f_out_msm);
fprintf('\nDispersion of the feature outside the MSM: %f',norerror);

if flag
  % representation des histogrammes de distribution de l'attribut
  % feature sur et en dehors de la MSM
  figure, subplot(2,2,1), imshow((msm~=C0).*feature,[]), colormap gray,
  title('Feature over the MSM');;
  subplot(2,2,2),    hist(f_on_msm,100);
  subplot(2,2,3), imshow((msm==C0).*feature,[]), colormap gray;
  title('Feature outside the MSM');
  subplot(2,2,4),    hist(f_out_msm,100);
  drawnow;    
end;

