%% SOURCE_ANGLE - Compute the orientation of a vectorial field.
%
%% Description
% Compute the orientation of a vectorial field. Used for the representation 
% of the source.
%
%% Syntax
%     source_vec  = source_angle( gx, gy )
% 
%% See also
% Related:     
% represent_source
% source

%% Function implementation
function source_vec = source_angle(gx, gy )

theta = angle( gx, gy );

% C0=255; CP=0; CM=127;
% c1 = (theta<=pi/4.);
% c2 = (theta>=7.*pi/4.);
% c3 = (theta<=5.*pi/4.);
% c4 = (theta<=3.*pi/4.);
% c5 = (~c1&~c2&~c3&~c4);
% source_vec =   (c1|c2 | c5) .* C0 +  (c3 | c5) .* CM ...
%     + (c4) .* CP;

source_vec=theta;

