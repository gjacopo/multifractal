%% BITS - Return the nearest upper powers of 2 of two values.
%
%% Syntax
%     [xeff, yeff] = bits( sx, sy)

%% Function implementation
function [xeff, yeff] = bits( sx, sy)

nxeff = floor(log(sx)/log(2.));
nyeff = floor(log(sy)/log(2.));

if sx~=2^nxeff nxeff=nxeff+1; end
if sy~=2^nyeff nyeff=nyeff+1; end

xeff = 2^nxeff;
yeff = 2^nyeff;
