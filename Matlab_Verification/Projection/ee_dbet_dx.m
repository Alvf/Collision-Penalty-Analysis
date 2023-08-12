function [dbdx] = ee_dbet_dx(numstate)
%EE_DBET_DX gets dbdx from collision projection frame
v0 = numstate(1:3);
v1 = numstate(4:6);
v2 = numstate(7:9);
v3 = numstate(10:12);

e0 = v1 - v0;
e1 = v3 - v2;
e2 = v2 - v0;

bet = (e1'*e1 * e2'*e0 - e2'*e1 * e0'*e1)/(e0'*e0 * e1'*e1 - e0'*e1 * e0'*e1);

e1n2 = e1'*e1;
e0n2 = e0'*e0;
e0d1 = e0'*e1; 
e0d2 = e0'*e2;
e1d2 = e1'*e2;
e0x12 = cross(e0,e1)'*cross(e0,e1);

dbdv0 = 2*bet*(e1n2*e0 - e0d1*e1) + e1*(e0d1 + e1d2) - e1n2*(e2 + e0);
dbdv1 = 2*bet*(e0d1*e1 - e1n2*e0) - e1*e1d2 + e1n2*e2;
dbdv2 = 2*bet*(e0n2*e1 - e0d1*e0) + e0*(e1n2 + e1d2) - 2*e0d2*e1 - e0d1*(e1 - e2);
dbdv3 = 2*bet*(e0d1*e0 - e0n2*e1) - e0*e1d2 + 2*e1*e0d2 - e0d1*e2;

dbdx = [dbdv0;dbdv1;dbdv2;dbdv3]/e0x12;

end

