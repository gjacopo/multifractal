%% CLIQUE - Extract a patch from an image.
%
%% Description
% Extraction of a cross-shaped clique:
%                        [iy-1][ix]  <- 3 
%    1 -> [iy][ix-1]      [iy][ix]  <- 0    [iy][ix+1]  <- 2
%                        [iy+1][ix]  <- 4
%
%% Syntax
%    parche = clique( img, ix, iy )

%% Function implementation
function parche = clique( img, ix, iy )

% [sx sy] = size(img);
% if(ix>1) iix=ix-1; else iix=sx; end;
% if(ix<sx) isx=ix+1; else isx=1;end;
% if(iy>1) iiy=iy-1; else iiy=sy;end;
% if(iy<sy) isy=iy+1; else isy=1;end;

parche = [img(ix,iy) img(ix-1,iy) img(ix+1,iy) ...
	  img(ix,iy-1) img(ix,iy+1)];
