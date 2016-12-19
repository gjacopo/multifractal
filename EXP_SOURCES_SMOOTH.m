dirpath='/nfs/data2/grazzini/Experience/fusion/Exp-fusion/';

% Nbre de bins de calcul de l'histogramme de representation 
% de la distribution des orientations: Nhisto
Nhisto=100;
% Nbre de bins a considerer pour cumuler les effectifs de pixels 
% dont l'orientation est dans [-s*pi, +s*pi]: nbbin
piprop=16;
s=1/piprop;
nbbin = floor(s / 2. * Nhisto);
s=2*s;
nbbin = [nbbin, floor(s / 2. * Nhisto)];

NBI=4;
NIM = cell(NBI,1);
NIM{1} = 'norma.ext.i20010624_S4.TC3';
NIM{2} = 'montana_landsat7';
NIM{3} = 'mississippi_landsat7';
NIM{4} = 'ext2.so.SPOT.2002_03_21_TC3';

NBJ=3;
EXT=cell(NBJ,1);
EXT{1}='_upm1.25';
EXT{2}='_upm1.5';
EXT{3}='_upm1.75';
% EXT{4}='';

orMM=cell(NBI*NBJ,1);

%% mississippi, upm 1.5
a=1
if a
for INDI=1:NBI 
  for INDJ=1:NBJ
    
    % nim='mississippi_landsat7';
    nim=NIM{INDI};
    cnim=[dirpath,nim,'.inr']
    img=load_inrimage(cnim,0,1);
    figure,    imagesc(img), colormap gray
 
    % ext='_upm1.5';
    ext=EXT{INDJ};
    ntest=[nim, ext];
    cntest=[dirpath,'contour.',ntest,'.inr']
    MSM=load_inrimage(cntest,0,1);
    
    [MM,orMMe,Mx,My,dummy] = ...
	source_frommsm( img, MSM, 1 );
    title(['Image ',ntest]);
    rec= reconstruction(dummy, Mx, My, flag);

    cnor=[dirpath,'or.',ntest,'.tif']
    imwrite(orMMe,jet,cnor,'tif');
    cnor=[dirpath,'or.',ntest,'.inr']
    % save_inrimage(cnor,orMMe,0,1,3);
    orMM{(INDI-1)*NBJ+INDJ} = orMMe;
    
    nf = [dirpath,'pctorient.',ntest,'.fi']
    [sx sy] = size(orMMe);
    [nn,xx] = hist(orMMe(:),Nhisto);
    nn = nn/ (sx*sy);
    figure, stem(xx,nn), 
    set(gca,'XTick',0:pi/2:2*pi)
    set(gca,'XTickLabel',{'0','pi/2','pi','3 pi/2','2 pi'})
    xlabel('0 \leq \Theta \leq 2\pi')
    ylabel('\delta_{\Theta}')
    xlim([0 2*pi])
    title(['Histo ',nf]);
    
    % Pct de pixels pour lesquels l'orientation est dans [-s*pi, s*pi]
    % et dans [-2*s*pi, 2*s*pi]
    nbinf=[];
    for i=1:2
      vinf = cumsum(nn(1:nbbin(i)));
      vsup = cumsum(nn(Nhisto-nbbin(i)+1:Nhisto));
      nbinf = [nbinf,vinf(nbbin(i)) + vsup(nbbin(i))]
    end;
    
    % fid = fopen(nf,'w');
    % fprintf(fid,'Pct pixels whom orientation is in [-pi/%d,pi/%d]: %f\n',...
    fprintf('Pct pixels whom orientation is in [-pi/%d,pi/%d]: %f\n',...
	    piprop, piprop, nbinf(1));
    % fprintf(fid,'Pct pixels whom orientation is in [-pi/%d,pi/%d]: %f\n',...
    fprintf('Pct pixels whom orientation is in [-pi/%d,pi/%d]: %f\n',...
	    piprop/2, piprop/2, nbinf(2));
    % fclose(fid);
    
  end;
end;

% save orMM;
end
return

NBI=2;
NIM{1} = 'norma.ext.MSAT6.1';
NIM{2} = 'lena256';
orMM_=cell(NBI,1);

for INDI=1:NBI 
  
  nim=NIM{INDI};
  cnim=[dirpath,nim,'.inr']
  img=load_inrimage(cnim,0,1);
  figure,    imagesc(img), colormap gray;
  
  ntest=nim;
  cntest=[dirpath,'contour.',ntest,'_upm1.25.inr']
  
  MSM=load_inrimage(cntest,0,1);
  
  figure,    imagesc(MSM), colormap gray
  
  [MM,orMMe,Mx,My,dummy] = ...
      source_frommsm( img, MSM, 1 );
  
  cnor=[dirpath,'or.',ntest,'.tif'];
  imwrite(orMMe,cnor,'tif');
  cnor=[dirpath,'or.',ntest,'.inr']
  % save_inrimage(cnor,orMMe,0,1,3);
  orMM_{INDI} = orMMe;
  
  [sx sy] = size(orMMe);
  [nn,xx] = hist(orMMe(:),Nhisto);
  nn = nn/ (sx*sy);
  figure, stem(xx,nn), 
  set(gca,'XTick',0:pi/2:2*pi)
  set(gca,'XTickLabel',{'0','pi/2','pi','3 pi/2','2 pi'})
  xlabel('0 \leq \Theta \leq 2\pi')
  ylabel('\delta_{\Theta}')
  xlim([0 2*pi])
  
  % Pct de pixels pour lesquels l'orientation est dans [-s*pi, s*pi]
  % et dans [-2*s*pi, 2*s*pi]
  nbinf=[];
  for i=1:2
    vinf = cumsum(nn(1:nbbin(i)));
    vsup = cumsum(nn(Nhisto-nbbin(i)+1:Nhisto));
    nbinf = [nbinf,vinf(nbbin(i)) + vsup(nbbin(i))]
  end;
  
  nf = [dirpath,'pctorient.',ntest,'.fi'];
  fid = fopen(nf,'w');
  fprintf(fid,'Pct pixels whom orientation is in [-pi/%d,pi/%d]: %f\n',...
	  piprop, piprop, nbinf(1));
  fprintf(fid,'Pct pixels whom orientation is in [-pi/%d,pi/%d]: %f\n',...
	  piprop/2, piprop/2, nbinf(1));
  fclose(fid);
end;

% save orMM_;


% Autre test: avec des bords calcules autrement,
% mais avec le meme noyau de reconstruction
nim='thres.cont.mshah.l15_norma.ext.i20010624_S4.TC3';
cnim=[dirpath,nim,'.inr'];
msm=load_inrimage(cnim,0,1);
MSM=(msm>60)*255;

nim='norma.ext.i20010624_S4.TC3';
cnim=[dirpath,nim,'.inr'];
img=load_inrimage(cnim,0,1);


[gx, gy] = derive_spectral( img );
[ax, ay] = derive_msm( MSM );

modo=0;
[gx, gy, MSM] = mask_gradient_msm( MSM, ax, ay, gx, gy, modo );
gx = propagation( gx, gy );
figure,imagesc(gx),colormap gray;

dummy =  dummy( MSM );

figure,imagesc(MSM),colormap gray;
figure,imagesc(dummy),colormap gray
