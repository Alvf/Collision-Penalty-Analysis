function [daldx] = ee_dal_dx(numstate)
%EE_DAL_DX gets dbdx from collision projection frame

v0 = numstate(1:3);
v1 = numstate(4:6);
v2 = numstate(7:9);
v3 = numstate(10:12);

e0 = v1 - v0;
e1 = v3 - v2;

al = e1'*e0/(e0'*e0);

e0n2 = e0'*e0;

dadv0 = 2*al*e0 - e1;
dadv1 = -dadv0;
dadv2 = -e0;
dadv3 = -dadv2;

daldx = [dadv0;dadv1;dadv2;dadv3]/e0n2;

end

