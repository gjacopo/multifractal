%% DUPLICATE - Duplicate the image edges to be periodic in both directions.
%
%% Syntax 
%      img = duplicate(img, width)

%% Function implementation
function img = duplicate(img, width)

[sx sy] = size(img);
img = [ img(sx-width+1:sx,:); img ; img(1:width,:)];
img = [img(:,sy-width+1:sy) img img(:,1:width)];
