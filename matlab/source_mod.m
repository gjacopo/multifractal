%% SOURCE_NORM - Compute the norm of a vectorial field.
%
%% Description
% Compute the norm of a vectorial field (eventually log-representation). 
% Used for the representation of the source.
%
%% Syntax
%      mod = source_norm( gx, gy, silog )
% 
%% See also
% Related:     
% represent_source
% source

%% Function implementation
function mod = source_norm( gx, gy, silog )

mod = sqrt( gx.^2 + gy.^2);
if silog
  mod = (mod>1e-17) .* mod + (mod<=1e-17) .* 10.^(-17.);
  mod = log(mod);
end;
