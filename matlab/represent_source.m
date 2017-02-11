%% REPRESENT_SOURCE - Represent the vectorial field of source.
%
%% Description
% Represent the vectorial field of source by its norm and its orientation.
%
%% Syntax
%     [source, source_vec] = represent_source( gx, gy, silog )
% 
%% See also
% Related:     
% source

%% Function implementation
function [source, source_vec] = represent_source( gx, gy, silog )


[sx sy] =size(gx);

source = source_norm( gx, gy, silog );
source_vec = source_angle(gx, gy);
