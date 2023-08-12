function [del] = del_dC_eedx(state)
%DEL_DC_EEDX Gets (dCee/dx)_Delta.
%   get vertices
v0 = state(1:3);
v1 = state(4:6);
v2 = state(7:9);
v3 = state(10:12);

e0 = v1-v0;
e1 = v3-v2;

z3 = zeros(3);

addpath('./Projection/')
dadx = ee_dal_dx(state);
dbetdx = ee_dbet_dx(state);
dgamdx = ee_dgam_dx(state);

dadv0 = dadx(1:3);
dadv1 = dadx(4:6);
dadv2 = dadx(7:9);
dadv3 = dadx(10:12);

dbdv0 = dbetdx(1:3);
dbdv1 = dbetdx(4:6);
dbdv2 = dbetdx(7:9);
dbdv3 = dbetdx(10:12);

dgdv0 = dgamdx(1:3);
dgdv1 = dgamdx(4:6);
dgdv2 = dgamdx(7:9);
dgdv3 = dgamdx(10:12);

del = [z3,z3,z3,z3;
       -e0*dadv0', -e0*dadv1', -e0*dadv2', -e0*dadv3';
       -e0*dbdv0'-e1*dgdv0', -e0*dbdv1'-e1*dgdv1', -e0*dbdv2'-e1*dgdv2', -e0*dbdv3' - e1*dgdv3'];

end

