%% QUANTILE_THRESHOLD - Compute the singularity exponent threshold value.
%
%% Description
% Compute the adaptative threshold giving the value of the singularity exponents
% of the pixels of the MSM, according to the desired density.
%
%% Syntax 
%     [umbral, n] = quantile_threshold( disp, upm_dens [, flag] )
%
%% Inputs
%    - disp : distribution of singularity exponents,
%    - upm_dens : desired density for the pixels of the MSM.
%
%% Outputs
%    - umbral : adaptative threshold for cuting the distribution
%      of exponents,
%    - n : true  density for the pixels of the MSM.
%
%% See also
% Related:    
% upm

%% Function implementation
function [umbral, n] = quantile_threshold( disp, upm_dens, flag ) 

if ~exist('flag')
  flag =0;
end;
  
[sx sy] = size(disp);
Nhisto=10000;

mind = min(min(disp));
maxd = max(max(disp));
wds = maxd - mind; 

if (wds>1e-30) 
  
  n = hist(disp(:),Nhisto) / (sx*sy);
  if flag
    % representation de l'histogramme de distribution des
    % exposants de singularite: on reduit le nombre de bins    
    
    nn = hist(disp(:),100);
    xx = mind:(wds/100.):maxd;
    xx=xx(2:101);
    figure, subplot(1,2,1), hist(disp(:),100), colormap (cool);
            subplot(1,2,2), stem(xx,nn), 
	    hold on, plot(xx,nn,'red'), hold off;
  end;
  
%  a = find(cumsum(n)>=(1.-upm_dens));
%  ip = a(1);
  
   ip=1; norma=0.;
   while (ip<Nhisto & norma<(1.-upm_dens))
     norma = norma + n(ip);
     ip=ip+1;
   end;
   umbral = maxd * ip /Nhisto;

else 
  umbral = mind;
end;
