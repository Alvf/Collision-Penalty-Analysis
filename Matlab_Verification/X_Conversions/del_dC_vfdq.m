function [delTerms] = del_dC_vfdq(state)
%DEL_DC_VFDQ gets (dCdq)_Delta from the state
v0 = state(1:3);
v1 = state(4:6);
v2 = state(7:9);
v3 = state(10:12);

e0 = v1 - v2;
e1 = v3 - v2;

zeros3 = zeros(3);

addpath('./Projection/');

alphdx = dal_dx(state);
betdx = dbet_dx(state);
gamdx = dgam_dx(state);

dadv1 = alphdx(4:6);
dadv2 = alphdx(7:9);
dadv3 = alphdx(10:12);

dbdv0 = betdx(1:3);
dbdv1 = betdx(4:6);
dbdv2 = betdx(7:9);
dbdv3 = betdx(10:12);

dgamdv0 = gamdx(1:3);
dgamdv1 = gamdx(4:6);
dgamdv2 = gamdx(7:9);
dgamdv3 = gamdx(10:12);

delTerms = [zeros3, zeros3, zeros3, zeros3;
            zeros3, -e0*dadv1', -e0*dadv2', -e0*dadv3';
            -e0*dbdv0'-e1*dgamdv0', -e0*dbdv1'-e1*dgamdv1', -e0*dbdv2'-e1*dgamdv2', -e0*dbdv3'-e1*dgamdv3'];

end

